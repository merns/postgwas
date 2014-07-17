readGWASdatasets <- function(
    gwas.datasets, 
    assert.positions.match = FALSE, 
    assert.unique.snps.per.dataset = FALSE,
    # the read.columns argument can specify several names of column headers for different file formats
    # each format has to contain the elements "SNP" and "P"
    read.columns = list( 
        PLINK = c(CHR = "CHR", BP = "BP", SNP = "SNP", P = "P"), 
        GEMMA = c(CHR = "chr", BP = "ps", SNP = "rs", P = "p_wald"), 
        FAsTLMM = c(CHR = "Chromosome", BP = "Position", SNP = "SNP", P = "Pvalue"),
        GenABEL = c(CHR = "Chromosome", BP = "Position", SNP = "snpnames", P = "P1df"), 
        custom = c(SNP = "SNP", P = "P")
    )
) {

  # check read.columns arg
  if(
      !is.list(read.columns) ||
      !all(sapply(read.columns, function(colset) all(c("SNP", "P") %in% names(colset)))) ||
      length(read.columns) < 1
  )
    stop("readGWASdatasets: Argument read.columns is malformed or does not contain elements 'SNP' and 'P'.\n")

  # check gwas.datasets arg (and make iterable, i.e. list or vector)
  if(is.data.frame(gwas.datasets) || class(gwas.datasets) == "scan.gwaa") {
    gwas.datasets <- list(dataset1 = gwas.datasets)
  } else if(is.list(gwas.datasets)) {
    if(is.null(names(gwas.datasets)))
      names(gwas.datasets) <- paste("dataset", 1:length(gwas.datasets), sep = "")
  } else if(is.character(gwas.datasets)) {
    names(gwas.datasets) <- gwas.datasets
  } else {
    stop("Malformed GWAS dataset argument, has to be a data frame, a character string (filename), an object of class 'scan.gwaa' or a list of these.")
  }
  if(any(duplicated(names(gwas.datasets))))
    stop("GWAS dataset names are not unique.")
  

  # loop over datasets and combine them
  message("Reading GWAS datasets ...")
  p.all <- NULL
  for( dsname in names(gwas.datasets) ) {
    
    ds <- gwas.datasets[[dsname]]
    
    # check for GenABEL object
    if(class(ds) == "scan.gwaa") {
      p.next <- cbind(snpnames = snpnames(ds), ds[, c("Chromosome", "Position", "P1df")])
    } else if(is.data.frame(ds) ) {
      p.next <- ds
    } else if(is.character(ds)) {
      p.next <- read.table(ds, header=TRUE)      
    } else {
      stop("Malformed GWAS dataset argument, has to be a data frame, a character string (filename), an object of class 'scan.gwaa' or a list of these.")
    }
    
    valid.formats <- which(sapply(read.columns, function(colset) all(colset %in% colnames(p.next))))
    if(length(valid.formats) > 0) {
      message(paste(names(read.columns)[valid.formats[1]], "format detected"))
      format <- read.columns[[valid.formats[1]]]
    } else {
      stop(paste(
              "Unknown format of gwas dataset '", dsname, "'. Recognized formats (column headers, arbitrary order):\n", 
              paste(colnames(read.columns), sapply(read.columns, paste, collapse = ", "), sep = ": ", collapse = "\n"),
              sep = ""
          ))
    }
    p.next <- p.next[, format]
    colnames(p.next) <- names(format)
    p.next$pheno = dsname
    
    p.next <- na.omit(p.next)
    if(assert.unique.snps.per.dataset && any(duplicated(p.next$SNP)))
      stop(paste("GWAS dataset '", dsname, "' contains duplicate SNPs (e.g. for PLINK files, check for multiple models in the file, e.g. DOM/REC/GENO... which is not allowed)."))
    
    if(!is.null(p.all) && (ncol(p.all) != ncol(p.next) || any(colnames(p.all) != colnames(p.next))))
      stop("Formats of two or more GWAS datasets are incompatible.\n")
    p.all <- rbind(p.all, p.next)
  }
  
  if(assert.positions.match && "BP" %in% names(format)) {
    p.all.test <- merge(p.all[, c("SNP", "BP")], p.all[, c("SNP", "BP")], by.x = "SNP", by.y = "SNP")
    if(any(p.all.test$BP.x != p.all.test$BP.y))
      stop(paste(
              "Base positions do not match between GWAS datasets for the following SNP (or the SNP occurs with different base positions in the same file): ",  
              as.vector(p.all.test[which(p.all.test$BP.x != p.all.test$BP.y)[1], "SNP"])
          ))
    rm(p.all.test)
  }
  
  return(p.all)
}

# returns data from a pedigree mapfile (linkage format, either 4 or 3 columns) without header. 
# Return value is a data frame with columns "CHR", "SNP", "BP"
readMapfile <- function(mapfile) {
  mapinfo <- read.table(mapfile, header = FALSE, row.names = NULL)
  if(ncol(mapinfo) == 4) mapinfo <- mapinfo[, c(1,2,4)]
  if(ncol(mapinfo) != 3) 
    stop("Mapfile format unknown")
  colnames(mapinfo) <- c("CHR", "SNP", "BP")
  return(mapinfo)
}
