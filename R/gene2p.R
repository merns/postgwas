# correction of association p by gene length (different number of SNPs per gene)
# most pathway tools use recalculation of the association statistic under a nulll model (phenotype label permutation)
# this requires genotype data and external tools and implies a dependency on certain study designs and test statitics
# such an implementation contradicts the aim of postgwas to be universally applicable and easy to set up
# alternative approches based on p-value lists only are rare or unsatisfactory
# thus, we use an alternative approach to derive a corrected assignment of SNP p-values to genes
# as outlined in the paper of Wang 2010, recombiantion hotspots, i.e. haplotype structure can be used to estimate the number of independet tests per gene, and this way a measure to correct the SNP-based pvalues. 
# We follow this suggestion by developing a

# The GATES version is recommended for case control studies
# it does not only consider LD, but allele freq and cohort size as well
# but as AF does not seem to impact correlation between SNPs (according to the original manuscript), it makes primarily sense for CC studies

# make sure arguments have to be sorted by increasing p-value
GATES <- function(ldmatrix, snps, p) {
  snpcount <- length(p)
  if(!all(dim(ldmatrix) == snpcount))
    stop("function GATES: Argument 'ldmatrix' is not rectangular or does not match the length of argument vector 'p'.\n")
  if(any(is.na(ldmatrix)))
    stop("function GATES: Argument 'ldmatrix' may not contain NA values.\n")
  if(snpcount < 1)
    stop("function GATES: No SNP provided.\n")
  if(snpcount == 1)
    # for a single SNP in a gene, there is a single independent P value:
    return(p)
  # this is the formula from the source code of GATES (class LDGeneBasedAssociationScan)
  ldmatrix <- 0.7723*ldmatrix^6 - 1.5659*ldmatrix^5 + 1.201*ldmatrix^4 - 0.2355*ldmatrix^3 + 0.2184*ldmatrix^2 + 0.6086*ldmatrix
  # remaining stuff according to Miao-Xin Li, 2011, AJHG 88:283-293
  eff.snpcount.fun <- function(ldmat) {
    ldmat <- as.matrix(ldmat)
    snpcount.local <- dim(ldmat)[1]
    if(snpcount.local <= 1)
      return(1)
    ev <- eigen(ldmat, only.values = TRUE)$values
    ev <- ev[ev > 1]
    snpcount.local - sum(ev - 1)
  }
  eff.snpcount.global <- eff.snpcount.fun(ldmatrix) 
  # modified Simes correction
  return(min(
    sapply(
      1:snpcount, 
      function(i) {
        (eff.snpcount.global * p[i]) / eff.snpcount.fun(ldmatrix[1:i, 1:i])
      }
    )
  ))
}

# calculates the number of effective independent SNPs  as SpD
# then uses the Sidak correction to adjust the p-values by that effective number of tests and returns the smallest
# this is slightly more conservative than the GATES procedure
SpD <- function(ldmatrix, snps, p) {
  snpcount <- length(p)
  if(!all(dim(ldmatrix) == snpcount))
    stop("function SpD: Argument 'ldmatrix' is not rectangular or does not match the length of argument vector 'p'.\n")
  if(any(is.na(ldmatrix)))
    stop("function SpD: Argument 'ldmatrix' may not contain NA values.\n")
  if(snpcount < 1)
    stop("function SpD: No SNP provided.\n")
  if(snpcount == 1)
    # for a single SNP in a gene, there is a single independent P value:
    return(p)
  ev <- eigen(ldmatrix, symmetric = TRUE, only.values = TRUE)$values
  # this is the formula suggested by Dale Nyholt 2004, AJHG 74:765-769
  eff.snpcount <- 1 + (snpcount -1) * (1 - var(ev) / snpcount)
  # return Sidak corrected p
  return(min(1 - (1 - p)^eff.snpcount))
}

# TODO gwas has to contain unique SNPs
gene2p <- function(
    gwas,
    gts.source, 
    method = GATES, 
    cores = 1
) {
  
  if(missing(gts.source))
    stop("gene2p: Argument 'gts.source' has to be specified")
  if(is.null(cores) || !is.numeric(cores) || cores <= 1) {
    cores = 1
  } else {
    require(parallel)
  }
  if(!is.function(method))
    stop("gene2p: Argument 'method' has to be a function object.\n")
  stopifnot(is.data.frame(gwas))
  stopifnot(c("SNP", "P") %in% colnames(gwas))
  if(!any(c("geneid", "genename") %in% colnames(gwas)))
    stop("gene2p: Argument 'gwas' neither contains a column 'geneid' nor 'genename'.\n")
  if(!is.function(method) || !all(c("ldmatrix", "snps", "p") %in% names(formals(method))))
    stop("gene2p: Argument 'method' has to be a function taking arguments ldmatrix, snps and p.\n")
  genecol <- which(colnames(gwas) %in% c("geneid", "genename"))[1]
  
  
  gwas <- vectorElements(gwas)
  # this is a list, named by genes and each element containing the vector of SNP identifiers of that gene
  snps.per.gene <- tapply(gwas$SNP, factor(gwas[, genecol]), function(x) as.vector(x), simplify = FALSE)
  snpcount.per.gene <- sapply(snps.per.gene, length)
  
  genos <- getGenotypes(snps = unlist(snps.per.gene), gts.source = gts.source, remove.homozygous = TRUE, toFile = "gene2p")
  
  # for each gene: 
  # make a list containing the following elements, all ordered by increasing p-value
  # ld: matrix of r square LD correlation values (columns / rows named, and ordered by SNP p)
  # snps: a vector of snps
  # p: vector of p
  calcLD <- function(snpset) {
    if(any(duplicated(snpset)))
      warning("gene2p: A SNP is annotated several times to the same gene. If the p-value is not identical for all occurences of a SNP, this will lead to an error.")
    # determine pairwise LD
    genotypes.avail <- which(genos@snpnames %in% snpset)
    if(length(genotypes.avail) > 1) {
      message(paste(" :", length(genotypes.avail)^2), appendLF = FALSE)
      rsqmatrix <- r2fast(data = genos, snpsubset = genotypes.avail)
      # rsqmatrix is symmetric matrix with r^2 values in the upper triangle and diag NA
      # eigen function uses the lower half and denies NA values: fill diag and lower triangle
      rsqmatrix[lower.tri(rsqmatrix)] <- t(rsqmatrix)[lower.tri(rsqmatrix)]
      diag(rsqmatrix) <- 1
      # return value: SNP order has to be by increasing p-value for all elements
      snps.p <- unique(gwas[gwas$SNP %in% colnames(rsqmatrix), c("SNP", "P")])
      snps.p <- vectorElements(snps.p[order(snps.p$P), ])
      rsqmatrix <- rsqmatrix[snps.p$SNP, snps.p$SNP]
      return(list(ld = rsqmatrix, snps = snps.p$SNP, p = snps.p$P))
    } else {
      if(length(genotypes.avail) <= 0) {
        return(NA)
      } else { # single genotype available
        snpnames.avail <- genos@snpnames[genotypes.avail]
        return(list(
                # do not use snps from snpset here, because that can be more than avail genotypes 
                ld = matrix(1, dimnames = list(snpnames.avail, snpnames.avail)),
                snps = snpnames.avail, 
                p = unique(gwas[gwas$SNP == snpnames.avail, "P"])
            ))
      }
    }
  }
  
  message("Calculating LD pairs for each gene", appendLF = FALSE)
  if(cores > 1) {
    ld.snp.p.per.gene <- mclapply(snps.per.gene, calcLD, mc.cores = cores)
  } else {
    ld.snp.p.per.gene <- lapply(snps.per.gene, calcLD)
  }
  message(" done.")
  
  # remove genes with NA data (e.g. no genotypes available)
  if(sum(is.na(ld.snp.p.per.gene)) > 0) {
    warning(paste("gene2p: ", sum(is.na(ld.snp.p.per.gene)), "genes were removed due to unavailable genotyp data (no SNPs with genotypes annotated to these genes).\n"))
    ld.snp.p.per.gene <- ld.snp.p.per.gene[!is.na(ld.snp.p.per.gene)]
  }
  # calculate the number of effective (i.e. independent) SNPs per gene, 
  # adjust the p-values and return a representative
  message("Calculating gene-wise p-values...")
  if(cores > 1) {
    gene.represent.p <- mclapply(ld.snp.p.per.gene, function(dat) method(dat$ld, dat$snps, dat$p), mc.cores = cores)
  } else {
    gene.represent.p <- lapply(ld.snp.p.per.gene, function(dat) method(dat$ld, dat$snps, dat$p))
  }
  res <- data.frame(names(gene.represent.p), gene.p = unlist(gene.represent.p))
  colnames(res)[1] <- colnames(gwas)[genecol]
  return(merge(gwas, res, all.x = TRUE))
}
