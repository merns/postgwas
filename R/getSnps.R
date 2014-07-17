getSnpsByWin <- function(chr, bp.start, bp.end, map = NULL, pop.id = 2) {
  if( is.null(chr) | is.null(bp.start) | is.null(bp.end) ) 
    stop("Position parameters (chr, bp.start, bp.end) not fully specified")

  if( is.null(map) & is.null(pop.id) ) 
    stop("Use either mapfile or Hapmap SNP retrieval (parameters 'map' or 'pop.id')")

  if( !is.null(map) ) {
    if(is.character(map)) map <- readMapfile(map)
    map <- vectorElements(map)
    return(vectorElements(map[map$CHR == chr & map$BP >= bp.start & map$BP <= bp.end, c("SNP", "CHR", "BP")]))
  } else {
    message("Retrieving HapMap SNPs")
    query.url <- url(
      paste("http://hapmap.ncbi.nlm.nih.gov/biomart/martservice?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20%20virtualSchemaName%20=%20%22rel27_NCBI_Build36%22%20formatter%20=%20%22TSV%22%20header%20=%20%220%22%20uniqueRows%20=%20%220%22%20count%20=%20%22%22%20datasetConfigVersion%20=%20%220.6%22%20%3E%20%20%20%20%20%20%20%20%3CDataset%20name%20=%20%22hm27_variation%22%20interface%20=%20%22default%22%20%3E%20%20%20%20%3CFilter%20name%20=%20%22stop%22%20value%20=%20%22", as.integer(bp.end), "%22/%3E%20%20%20%20%3CFilter%20name%20=%20%22chrom%22%20value%20=%20%22chr", chr, "%22/%3E%20%20%20%20%3CFilter%20name%20=%20%22pop_id%22%20value%20=%20%2", pop.id, "2%22/%3E%20%20%20%20%3CFilter%20name%20=%20%22start%22%20value%20=%20%22", as.integer(bp.start), "%22/%3E%20%20%20%20%3CAttribute%20name%20=%20%22marker1%22/%3E%20%20%20%20%3CAttribute%20name%20=%20%22chrom%22%20/%3E%%3CAttribute%20name%20=%20%22start%22%20/%3E%20%20%3C/Dataset%3E%3C/Query%3E", sep = "")
    )
    res <- read.table(query.url, col.names = c("SNP", "CHR", "BP"))
    # HapMap requires numerical chromosomes e.g. "3" for queries, but returns character chromosome names like "chr3" ... this does not make sense
    res$CHR <- sub("chr", "", res$CHR, fixed = TRUE)
    return(vectorElements(res))
  }

}


getSnpsByRS <- function(snps, map = NULL, pop.id = 2) {

  if( is.null(snps) || !is.character(snps) ) 
    stop("Parameter snps not specified or in wrong format\n")

  if( is.null(map) & is.null(pop.id) ) 
    stop("Use either mapfile or Hapmap SNP retrieval (parameters 'map' or 'pop.id')\n")

  if( !is.null(map) ) {

    if(is.character(map)) map <- readMapfile(map)
    snps.mapped <- map[map$SNP %in% snps, c("SNP", "CHR", "BP")]
    
  } else {

    # map snps to correct base positions
		snps.mapped <- data.frame(SNP = NULL, CHR = NULL, BP = NULL)
		snps.blocks <- makeBlocks(snps, block.size = 250)
		for( snps.block in snps.blocks ) {
			message("Retrieving HapMap SNPs")
			query.url <- url(paste("http://hapmap.ncbi.nlm.nih.gov/biomart/martservice?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20%20virtualSchemaName%20=%20%22rel27_NCBI_Build36%22%20formatter%20=%20%22TSV%22%20header%20=%20%220%22%20uniqueRows%20=%20%220%22%20count%20=%20%22%22%20datasetConfigVersion%20=%20%220.6%22%20%3E%20%20%20%20%20%20%20%20%3CDataset%20name%20=%20%22hm27_variation%22%20interface%20=%20%22default%22%20%3E%20%20%20%20%3CFilter%20name%20=%20%22marker_name%22%20value%20=%20%22", paste(snps.block, collapse = ","), "%22/%3E%3CFilter%20name%20=%20%22pop_id%22%20value%20=%20%22", pop.id, "%22/%3E%20%20%20%20%3CAttribute%20name%20=%20%22marker1%22/%3E%20%20%20%20%3CAttribute%20name%20=%20%22chrom%22%20/%3E%%3CAttribute%20name%20=%20%22start%22%20/%3E%20%20%3C/Dataset%3E%3C/Query%3E", sep = ""))
			snps.mapped <- rbind(snps.mapped, read.table(query.url, col.names = c("SNP", "CHR", "BP")))
		}
    # HapMap requires numerical chromosomes e.g. "3" for queries, but returns character chromosome names like "chr3" ... this does not make sense
    snps.mapped$CHR <- sub("chr", "", snps.mapped$CHR, fixed = TRUE)
		
  }
  
  if(!all(snps %in% snps.mapped$SNP))
    stop(paste("Cannot retrieve position information from the given source (HapMap or genotype file) for the following SNPs (remove from source data?):", paste(snps[!(snps %in% snps.mapped$SNP)], collapse = ",")))
  
  return(vectorElements(snps.mapped))

}
