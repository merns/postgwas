ldGenes <- function(snps, genes, gts.source = 2, ld.win = 1000000, cores = 1) {
  
  if(!is.numeric(snps$BP) | !is.numeric(genes$start) | !is.numeric(genes$end))
    stop("Base positions have to be numeric")

  gts.source.parsed <- parseGtsSourceArg(gts.source)
  
  # pre-read map data
  if(gts.source.parsed$mode == "gwaa-files" || gts.source.parsed$mode == "gwaa-object") {
    if(gts.source.parsed$mode == "gwaa-files")
      capture.output(gwaa <- load.gwaa.data(phenofile = gts.source.parsed$phe.fn, genofile = gts.source.parsed$gwaa.fn))
    else
      gwaa <- gts.source.parsed$snp.data.obj
    rm(gts.source.parsed)
    map.dat <- data.frame(CHR = chromosome(gwaa), SNP = snpnames(gwaa), BP = map(gwaa))
    rm(gwaa)
  } else if(gts.source.parsed$mode == "linkage") {
    map.dat <- readMapfile(gts.source.parsed$map.fn)
  } else {
    map.dat <- NULL
  }
  
  # get positions of SNPs according to genotype source
  # snps$BP matches genes positions, snps$BP.mapped matches genotype positions
  snps.mapped <- getSnpsByRS(as.vector(snps$SNP), map = map.dat)
  snps <- vectorElements(merge(snps, snps.mapped, by = "SNP", suffixes = c("", ".mapped")))
  genes <- genes[order(genes$start, na.last = NA), ]
    
  # get genes for each SNP (within window + nearest)
  querysnp.genes <- lapply(
    1:nrow(snps), 
    function(idx) {
      bp <- snps[idx, "BP"]
      genes.chr <- genes[genes$chrname == snps[idx, "CHR"], ]
      genes.chr$CHR.mapped <- snps[idx, "CHR.mapped"] 
      nearest.bp <- min(abs(genes.chr$start - bp), abs(genes.chr$end - bp))
      return(unique(rbind(
        # window genes
        genes.chr[genes.chr$end >= (bp - ld.win/2) & genes.chr$start <= (bp +ld.win/2), ], 
        # nearest gene
        genes.chr[abs(genes.chr$start - bp) == nearest.bp | abs(genes.chr$end - bp) == nearest.bp, ] 
      )))
     }
  )
  names(querysnp.genes) <- snps$SNP
  
  
  # define function to collect snps in genes
  # returns a nested list: querysnps -> genes -> snps in gene  
  # e.g returnedlist$querysnp$genename contains rs-ids of all snps in that gene (max. 100)
  collectSnpsPerGene <- function(querysnp.idx) { 
    querysnp <- as.vector(snps$SNP)[querysnp.idx] 
    shift <- as.vector(snps$BP - snps$BP.mapped)[querysnp.idx]
    genes <- querysnp.genes[[querysnp]]

    # get snps within window (when map.dat is NULL, uses HapMap retrieval)
    ldsnps.all <- getSnpsByWin(
        genes$CHR.mapped[1], 
        min(genes$start - shift), 
        max(genes$end - shift), 
        map = map.dat
    )
    # return only snps within genes
    return(mapply(
            function(id, start, end) {
              ldsnps.gene <- ldsnps.all[as.vector(ldsnps.all$BP) >= start & as.vector(ldsnps.all$BP) <= end, ]
              ldsnps.gene <- as.vector(ldsnps.gene[order(ldsnps.gene$BP), "SNP"])
              return(pruneVec(ldsnps.gene, 100))
            }, 
            as.vector(genes$genename),  # dummy, only used for mapply to return a named list
            genes$start - shift, 
            genes$end - shift,
            SIMPLIFY = FALSE
        ))
  }
  
  # run function to collect snps in genes
  if(cores > 1) {
    require(parallel)
    querysnp.genes.snps <- mclapply(1:length(snps$SNP), collectSnpsPerGene, mc.cores = cores)
  } else {
    querysnp.genes.snps <- lapply(1:length(snps$SNP), collectSnpsPerGene)
  }
  names(querysnp.genes.snps) <- snps$SNP
  
  
  message("Loading Genotypes...")
  snps.all <- unique(c(
    unlist(querysnp.genes.snps),    # snps from genes
    names(querysnp.genes.snps)      # always include the query SNP itself
  ))
  genos <- getGenotypes(snps = snps.all, gts.source = gts.source, remove.homozygous = TRUE, toFile = "snp2gene")

  # define function to calc LD for all query SNPs for all assigned genes
  calcLD <- function(querysnp) {
    querysnp.avail <- querysnp %in% genos@snpnames
    genes.snps <- querysnp.genes.snps[[querysnp]]
    ld <- sapply(
        names(genes.snps), 
        function(gene) {
          snps <- genes.snps[[gene]]
          querysnp.incl <- querysnp %in% snps
          # if query SNP among SNPs in gene, remove that (r2fast cannot calculate LD between identical SNPs)
          if(querysnp.incl)
            snps <- snps[snps != querysnp]
          snps.genos.avail <- genos@snpnames[genos@snpnames %in% snps]
          if(!querysnp.avail || length(snps.genos.avail) < 1)
            return(c(ld.max = 0, ld.mean = 0, ld.sdev = 0))
          message(paste("Calculating LD:", querysnp, "<->", gene, "(", length(snps.genos.avail), "snps )"))
          ld <- as.vector(r2fast(data = genos, snpsubset = snps.genos.avail, cross.snpsubset = querysnp)$r^2)
          if(querysnp.incl)
            c(ld, 1)
          return(c(
                  ld.max = round(max(ld), digits = 4),
                  ld.mean = round(mean(ld), digits = 4),
                  ld.sdev = round(sd(ld), digits = 4)
          ))
        }
    )
    # join the source data frame 'snps' with the resulting LD information (querysnp.genes contains data from 'snps' argument plus tested genes) 
    ldgenes.res <- as.data.frame(cbind(querysnp.genes[[querysnp]], t(ld)))
    ldgenes.res$SNP <- querysnp
    return(ldgenes.res)
  }
  
  # run LD calculation
  if(cores > 1) {
    require(parallel)
    res <- mclapply(snps$SNP, calcLD, mc.cores = cores)
  } else {
    res <- lapply(snps$SNP, calcLD)
  }

  return(list2df(res))
  
}

