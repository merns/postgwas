removeNeighborSnps <- function(
                        snps, 
                        maxdist = 50000, 
                        biomart.config = biomartConfigs$hsapiens,
                        use.buffer = FALSE
) {
  
  if(is.null(biomart.config) || !is.list(biomart.config))
    stop("Argument 'biomart.config' has to be a list")                        
                        
  if(nrow(snps) < 1)
    return(snps)
  
  if(!is.null(snps$BP) & !is.null(snps$CHR)) {
    # dont need further information
  } else if(!is.null(snps$SNP)) {
    message("Fetching SNP information from biomart...")
    snps.bm <- bm.snps(filter.val = snps$SNP, config = biomart.config, use.buffer = use.buffer)
    names(snps.bm)[ names(snps.bm) == "chrname" ] <- "CHR"
    names(snps.bm)[ names(snps.bm) == "chrom_start" ] <- "BP"
    # merge will sort the data frame, SNPs that do not match with ensembl are removed
    snps <- merge( snps, snps.bm, by.x = "SNP", by.y = "refsnp_id", suffixes = c(".original","") )
  } else {
    stop("Parameter 'snps' has to contain either column 'SNP' or 'BP' and 'CHR'")
  }
  
  snps <- snps[!is.na(snps$BP) & !is.na(snps$CHR), ]
  
  if(is.null(snps$P)) {
    # we require SNPs not to have identical BPs and CHR (would both be removed)
    snps <- snps[!duplicated(snps[, c("CHR", "BP")]), ]
    # for each SNP, remove when there are SNPs within maxdist AND having larger index
    return(snps[
              sapply(
                1:nrow(snps),
                function(idx) {
                  # snps in range
                  cand.idx <- which(
                      as.vector(snps$CHR) == as.vector(snps[idx, "CHR"]) & 
                      abs(as.numeric(snps$BP) - as.numeric(snps[idx, "BP"])) <= maxdist
                    )
                  if(length(cand.idx) < 1) 
                    return(TRUE)
                  if(any(max(cand.idx) > which(as.vector(snps$CHR) == as.vector(snps[idx, "CHR"]) & as.numeric(snps$BP) == as.numeric(snps[idx, "BP"]))))
                    return(FALSE)
                  else
                    return(TRUE)
                }
                )
              , ])
  } else {
    # we do that recursively:
    # take index with lowest p
    # invalidate all indices within window range
    # proceed with next lowest p that is valid and so on
    
    # recursive helper function 
    winsel<- function(bp, p) {
      winsel.intern <- function(bp, idx.valid) {
        if(any(idx.valid)) {
          idx <- min(which(idx.valid))
          # invalidate all SNPs within the window of idx 
          idx.valid[abs(bp - bp[idx]) < maxdist] <- FALSE
          return(c(idx, winsel.intern(bp, idx.valid)))
        } else {
          return()
        }
      }
      
      if(!is.vector(bp) || !is.vector(p) || length(bp) !=  length(p) || length(bp) < 1)
        stop("Arguments to the winsel function are not vectors or do not have equal length > 0")
      if(any(bp != bp[order(p)]))
        stop("The bp vector has to be ordered by p descendingly.")
      bp <- as.numeric(bp)
      sel.idx <- winsel.intern(bp, rep(TRUE, length(bp)))
    }
    
    # process chromosomes one by one
    do.call(
      rbind, 
      by(
        snps, 
        snps$CHR, 
        function(bychr) {
          bychr <- bychr[order(bychr$P), ]
          bychr[winsel(bychr$BP, bychr$P), ]
        }
      )
    )

  }
}
