data.regionalplot.ld <-function(regions.df, ld.options, cores) {
  
  message("Extracting SNPs for LD computation: ")
  # for each region, remove snps so that we have at most max.snps.per.window per region (evenly selected)
  regions.df <- regions.df[order(regions.df$BP.mapped), ]
  # a vector of rs ids
  ldsnps <- unique(unlist(tapply(
              regions.df$SNP, 
              factor(regions.df$packet.snp), 
              pruneVec, 
              ld.options$max.snps.per.window
          )))
  
  # this will be a snp.data object of GenABEL
  genos <- getGenotypes(
      snps = ldsnps, 
      gts.source = ld.options$gts.source, 
      remove.homozygous = TRUE,
      toFile = "regionalplot"
  )
  
  # function to calculate LD matrices for each region, to be called later
  calcLD <- function(packet.snp) {
    region.snps <- as.vector(regions.df$SNP[regions.df$packet.snp == packet.snp])
    # select genos that belong to snps from that region
    genotypes.avail <- which(genos@snpnames %in% unique(region.snps))
    if(length(genotypes.avail) < 2) return(NULL)
    message("Precalculating ", length(genotypes.avail)^2, " LD pairs for current region")
    rsq <- r2fast(data = genos, snpsubset = genotypes.avail)
    rsq[lower.tri(rsq)] <- NA
    return(rsq)
  }
  
  if(cores > 1) 
    rsq <- mclapply(unique(regions.df$packet.snp), calcLD, mc.cores = cores)
  else
    rsq <- lapply(unique(regions.df$packet.snp), calcLD)

  names(rsq) <- unique(regions.df$packet.snp)
  return(rsq)
}