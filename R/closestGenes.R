closestGenes <- function(snps, genes, mode = "cover", level = 0) {
  # For all snps, gets the 'level'-next genes (up- downstream or within)
  # 
  # Args:
  #  snps:         A data frame having columns 'SNP', 'CHR' and BP for reference SNP IDs, chromosome and base position (numeric)
  #  genes:        A data frame having numeric columns "chrnames", "start" and "end" (may have additional columns like name or ID).
  #                May contain additional columns (e.g. gene symbol) which are preserved in the return value.
  #  mode:         can be "cover", "up" or "down" to retrieve the next gene in the named direction. 
  #                "up" means smaller (gene end) base position in comparison to the SNP. 
  #                "down" means larger gene start.
  #                "cover" means smaller gene start and larger gene end.
  #  level:        0 refers to the closest gene in the direction given by mode, everything else the 'level' next gene in that direction.
  #                For covering genes, proximity is determined by the "start" gene column in comparison to the SNP.
  #
  # Value:
  #  A data frame with rows according to the query SNPs, containing further columns with information for the next gene (taken from 'genes' parameter)
  #  For SNPs that lack a next gene at the requested level for any reason, row values are set to NA.

  if(!is.numeric(snps$BP) | !is.numeric(genes$start) | !is.numeric(genes$end))
    stop("Base positions have to be numeric")

  # use mapping to generate numeric chromosomes (order and number does not matter)
  chrnames <- unique(c(genes$chrname, snps$CHR))
  chr.map  <- data.frame(chrname = chrnames, chromo = 1:length(chrnames))
  genes    <- vectorElements(merge(genes, chr.map, all.x = T))
  snps     <- vectorElements(merge(snps, chr.map, by.x = "CHR", by.y = "chrname", all.x = T))
  
  # make base positions single sequence for all chromosomes by setting pos = pos + chrom * shift
  base.shift <- 10^ceiling(log10(max(c(genes$end, snps$BP))))  # this is the number of digits needed for the largest number of bases - further digits can be used to make chromos sequential
  genes$start <- genes$chromo * base.shift + genes$start
  genes$end   <- genes$chromo * base.shift + genes$end
  snps$BP <- snps$chromo * base.shift + snps$BP

  if( mode == "up" ) {
    genes <- genes[order(genes$end, na.last = NA), ]             # sort
    genes.idx <- findInterval(snps$BP, genes$end) -level         # for all snps, gives indices in genes.sort where geneend smaller than bp of snp
    genes.idx[genes.idx <= 0 | genes.idx > nrow(genes)] <- NA    # out of bounds indices to NA
    res <- genes[genes.idx, ]
	
  } else if( mode == "down" ) {
    genes <- genes[order(genes$start, na.last = NA), ]           # sort
    genes.idx <- findInterval(snps$BP, genes$start) +level +1    # for all snps, gives indices in genes.sort where genestart larger than bp of snp
    genes.idx[genes.idx <= 0 | genes.idx > nrow(genes)] <- NA    # out of bounds indices to NA
    res <- genes[genes.idx, ]
	
  } else if( mode == "cover" ) {
    # cover gene does not necessarily need to be the next up or down gene. We have to check all genes here
    res <- list2df(lapply(
      snps$BP, 
      function(snp.bp) {
        cover.cands <- genes[genes$start < snp.bp & genes$end > snp.bp, ]
        if(nrow(cover.cands) > 1) {
          # order by average dist from SNP
          cover.cands <- cover.cands[order(cover.cands$start + cover.cands$end - 2 * snp.bp), ]
          return(cover.cands[level +1, ])
        } else {
          return(cover.cands[level +1, ]) # returns NAs when nrow(cover.cands) == 0
        }
      }
    ))
    
    # this is the fast but incorrect version for cover
#    genes <- genes[order(genes$start, na.last = NA), ]           # sort
#    genes.idx <- findInterval(snps$BP, genes$start) -level       # for all snps, gives indices in genes.sort where genestart smaller than bp of snp
#    genes.idx[genes.idx <= 0 | genes.idx >= nrow(genes)] <- NA   # out of bounds indices to NA
#    res <- genes[genes.idx, ]                                    # genes that start before snp
#    res[na.set(res$end <= snps$BP), ] <- NA                      # and end behind snp
  } 
   

  # join with snps data (both have same nrow after findInterval)
  res <- data.frame(res, snps, direction = mode)
  # de-sequentialize base postiions
  res$start <- res$start - res$chromo * base.shift
  res$end   <- res$end - res$chromo * base.shift
  res$BP <- res$BP - res$chromo * base.shift
  # remove annotations beyond chromosomes that arose through sequentialization
  # i.e. check that genes and SNPs chromosomes match
  res[na.set(res$chrname != res$CHR), ] <- NA

  return(res)
}
