  # This funtion can be used to complete the mapping of base positions between different assemblies.
  # When for a list of SNPs the original position is given and the mapped position is incomplete, 
  # the next position of the original assemly will be used to calculate the offset between the assemblies.
  # This offest is then used to impute the new unknown position of the neighbor SNP.
  # the data frame will be sorted by its BP column during this process. 
  # snps: a data frame containing the columns "BP"and "BP.mapped". The BP.mapped column will be imputed at the positions where it is NA. Attention: At least one cell of BP.mapped has to contain a valid position!
  # remove.unmappable: when FALSE, stops when snps cannot be mapped, otherwise returns an empty dataframe (all unmappable snps removed)
imputeRemappedPositions <- function(snps, remove.unmappable = FALSE) {
  
  if(is.null(snps) || nrow(snps) < 1)
    return(snps)
  
  if(!all(c("BP", "BP.mapped") %in% colnames(snps)))
    stop("Invalid argument: Data frame must contain the columns BP and BP.mapped")
  
  snps <- snps[order(as.numeric(snps$BP)), ]
  
  snps.idx.fill <- which(is.na(snps$BP.mapped))
  snps.idx.ok <- which(!is.na(snps$BP.mapped))
  
  # assert there is at least one NA and one non-NA index
  if(length(snps.idx.fill) < 1) 
    return(snps)
  if(length(snps.idx.ok) < 1) {
    if(remove.unmappable) {
      warning(paste("Removing SNPs with unknown base positions that cannot be imputed (no reference position on same chromosome): ", paste(sapply(snps, paste, collapse = ","), collapse = " - ")))
      return(na.omit(snps))
    }
    stop("Cannot impute unknown SNP base positions - no reference position available")
    return(na.omit(snps))
  }
  
  message(paste("Position of", length(snps.idx.fill), "SNPs are unkown in biomart or source data. Trying to estimate by neighbor offset."))
  
  # this code finds for each NA-index the next index that is not NA
  snps.idx.next.ok <- sapply(
    snps.idx.fill, 
    function(fill.idx) snps.idx.ok[which.min(abs(snps.idx.ok - fill.idx))]
  )
  
  # populate the NA cells
  snps[snps.idx.fill, "BP.mapped"] <- snps[snps.idx.fill, "BP"] - (snps[snps.idx.next.ok, "BP"] - snps[snps.idx.next.ok, "BP.mapped"])
  return(snps)
}
