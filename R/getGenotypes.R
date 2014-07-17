getGenotypes <- function(
  snps,
  gts.source = 2, 
  remove.homozygous = FALSE,
  return.characters = FALSE,
  return.genotypes = FALSE,
  toFile = "getGenotypes"
){
  
  ###################################
  ######### PARAMETER CHECK #########
  ###################################
  
  if(missing(snps))
    stop("getGenotypes: Parameter 'snps' has to be specified")
  
  if(is.factor(snps))
    snps <- as.vector(snps)
  if(!is.vector(snps))
    stop("getGenotypes: Parameter 'snps' has to be a character vector\n")
  
  if(is.null(snps) | length(snps) < 1) {
    warning("getGenotypes: No SNPs specified to retrieve genotypes for\n")
    return(NULL)
  }
  
  gts.source <- parseGtsSourceArg(gts.source)
  
  # check toFile argument
  if(is.null(toFile)) {
    if(gts.source$mode %in% c("linkage", "hapmap")) 
      stop("getGenotypes: Argument 'toFile' can only be NULL when 'gts.source' is a GenABEL gwaa data file or object.")    
    else 
      toFile = FALSE    
  } else {
    if(is.vector(toFile) && is.character(toFile) && length(toFile) == 1) {
      base.out.fn <- toFile
      ped.out.fn <- paste(toFile, "ped", sep = ".")
      map.out.fn <- paste(toFile, "map", sep = ".")
      gwaa.out.fn <- paste(toFile, "gwaa", sep = ".")
      phe.out.fn <- paste(toFile, "phe", sep = ".")
      toFile = TRUE
    } else {
      stop("getGenotypes: Malformed 'toFile' argument.")
    }    
  }
  
  
  
  ####################################
  ####### MODE GWAA (GenABEL)  #######
  ####################################
  if( gts.source$mode == "gwaa-files" || gts.source$mode == "gwaa-object" ) {
    
    if(gts.source$mode == "gwaa-object") 
      gts <- gts.source$snp.data.obj
    else
      capture.output(gts <- load.gwaa.data(phenofile = gts.source$phe.fn, genofile = gts.source$gwaa.fn)@gtdata)
    # gts is now an object of class snp.data
    
    snp.colmask <- gts@snpnames %in% snps
    
    if(!any(snp.colmask)) {
      warning("getGenotypes: None of the requested SNPs are contained in the genotype dataset provided.")
      return(NULL)
    }
    
    # extract SNPs
    gts <- gts[, snp.colmask]
    
    if(remove.homozygous) {
      maf.snps <- summary(gts)$Q.2
      if(any(maf.snps == 0))      # conditional: avoid subsetting with zero selection (throws error)
        gts <- gts[, maf.snps != 0]
    }

    # seems to me that export.plink does not do what I expected it to do ... better not use it
#    if(toFile) {
#      message(paste("Writing genotypes to ", gwaa.out.fn, " / ", phe.out.fn, " / ", ped.out.fn, " / ", map.out.fn, " ...", sep = ""))
#      capture.output(export.plink(data = gwaa, filebasename = base.out.fn, phenotypes = NULL))
#      capture.output(save.gwaa.data(data = gwaa, phenofile = phe.out.fn, genofile = gwaa.out.fn))
#    }
    
  } else { 
    ###################################
    ########### OTHER MODES ###########
    ###################################
    
    gts <- NULL  # this will be a matrix of character genotypes after retrieval
    mapinfo <- NULL # this will be a data frame with columns "CHR", "SNP", "BP"
    
    if( gts.source$mode == "linkage" ) {
      ###################################
      ####### PED / MAP RETRIEVAL #######
      ###################################
      message("Reading LINKAGE genotype files ...")
      mapinfo <- readMapfile(gts.source$map.fn)
      snp.colmask <- mapinfo$SNP %in% snps
      
      if(!any(snp.colmask)) {
        warning("getGenotypes: No SNPs were extracted from ped file... is any of the requested SNPs contained in the mapfile?\n")
        return(NULL)      
      }
      
      # read ped file in chunks, extract SNPs of interest
      pedfile <- file(gts.source$ped.fn, open = "r")
      
      while(length(ped <- readLines(pedfile, n = 35)) > 0) {
        line.elem.count <- length(unlist(strsplit(ped[1], split = "\\s+")))
        ped <- unlist(strsplit(ped, split = "\\s+"))  # yields a single vector, each line from pedfile split and concatenated with following line
        # single row will also be a matrix
        ped <- matrix(ped, ncol = line.elem.count, byrow = TRUE) # ped file as matrix, each field in file one matrix element
        # remove leading non-genotype columns (always select > 1 column, even single SNP has 2 cols for each allele, thus cannot be coerced to vector)
        ped <- ped[, (line.elem.count - 2*length(snp.colmask) +1):line.elem.count]
        # select target snps and combine the two alleles for each SNP
        ped <- paste(ped[, (which(snp.colmask) *2) -1], ped[, which(snp.colmask) *2], sep = "") # returns a long vector (concatenated by columns (snps))
        gts <- rbind(gts, matrix(data = ped, ncol = sum(snp.colmask)))
        rm(ped, line.elem.count)
      }
      close(pedfile)
      
      message("Formatting genotype data...")
      # format for 'genotype' class (recode numeric nucleotides to character)
      # gts is here automatically coerced to a vector
      gts <- gsub("1", "A", gts, fixed = TRUE)
      gts <- gsub("2", "C", gts, fixed = TRUE)
      gts <- gsub("3", "G", gts, fixed = TRUE)
      gts <- gsub("4", "T", gts, fixed = TRUE)
      if(is.null(gts) || length(gts) < 1) {
        warning("getGenotypes: No SNPs were extracted from ped file... is any of the requested SNPs contained?\n")
        return(NULL)
      }
      # convert to matrix, with SNPs listed in columns
      gts <- matrix(data = gts, ncol = sum(snp.colmask))
      colnames(gts) <- mapinfo$SNP[snp.colmask]
      
    } else if( gts.source$mode == "hapmap" ) {
      ###################################
      ####### HAPMART RETRIEVAL #########
      ################################### 
      # hapmap queries do not work with r package 'biomart' because the virtualSchemaName is not properly set (always 'default' instead of e.g. 'rel27_NCBI_Build36')
      # construct the RESTful XML query address for martservice
      # number of SNPs might be large, divide into 250-element blocks
      snps.blocks <- makeBlocks(snps, block.size = 250)
      for( snps.block in snps.blocks ) {
        message("Retrieving genotypes from HapMart...")
        query.url <- url(paste("http://hapmap.ncbi.nlm.nih.gov/biomart/martservice?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20%20virtualSchemaName%20=%20%22rel27_NCBI_Build36%22%20formatter%20=%20%22TSV%22%20header%20=%20%220%22%20uniqueRows%20=%20%220%22%20count%20=%20%22%22%20datasetConfigVersion%20=%20%220.6%22%20%3E%3CDataset%20name%20=%20%22hm27_variation%22%20interface%20=%20%22default%22%20%3E%3CFilter%20name%20=%20%22marker_name%22%20value%20=%20%22", paste(snps.block, collapse = ","), "%22/%3E%3CFilter%20name%20=%20%22pop_id%22%20value%20=%20%22", gts.source$pop.id, "%22/%3E%22/%3E%3CAttribute%20name%20=%20%22marker1%22%20/%3E%3CAttribute%20name%20=%20%22genotype%22%20/%3E%3C/Dataset%3E%3C/Query%3E", sep = ""))
        gts <- c(gts, readLines(query.url))
        close(query.url)
        if(!is.null(toFile)) {
          # download base positions for mapfile
          query.url <- url(paste("http://hapmap.ncbi.nlm.nih.gov/biomart/martservice?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20%20virtualSchemaName%20=%20%22rel27_NCBI_Build36%22%20formatter%20=%20%22TSV%22%20header%20=%20%220%22%20uniqueRows%20=%20%220%22%20count%20=%20%22%22%20datasetConfigVersion%20=%20%220.6%22%20%3E%3CDataset%20name%20=%20%22hm27_variation%22%20interface%20=%20%22default%22%20%3E%3CFilter%20name%20=%20%22marker_name%22%20value%20=%20%22", paste(snps.block, collapse = ","), "%22/%3E%3CFilter%20name%20=%20%22pop_id%22%20value%20=%20%22", gts.source$pop.id, "%22/%3E%22/%3E%3CAttribute%20name%20=%20%22chrom%22%20/%3E%3CAttribute%20name%20=%20%22marker1%22%20/%3E%%3CAttribute%20name%20=%20%22start%22%20/%3E%3C/Dataset%3E%3C/Query%3E", sep = ""))
          mapinfo <- rbind(mapinfo, read.table(query.url, col.names = c("CHR", "SNP", "BP")))
          # HapMap requires numerical chromosomes e.g. "3" for queries, but returns character chromosome names like "chr3" ... this does not make sense
          mapinfo$CHR <- sub("chr", "", mapinfo$CHR, fixed = TRUE)
        }
      }
      if(length(gts) < 1) {
        warning("getGenotypes: No SNPs were returned from HapMart... is human SNP data being used?\n")
        return(NULL)
      }
      message("Formatting genotype data...")
      # extract rs numbers and genotypes as strings
      gts <- strsplit(gts, split = "\t", fixed = TRUE) 
      gts <- unlist(gts)
      gts.rsnumbers <- gts[seq(1, length(gts), by = 2)]
      gts <- gts[seq(2, length(gts), by = 2)]
      # convert genotype strings to vectors (containing single genotypes per individ)
      gts <- sapply(gts, strsplit, split = " ", fixed = TRUE)
      # make genotype vectors uniform length (appending NAs by selecting over-length subset)
      gts.maxcount <- max(sapply(gts, length))
      gts <- array(gts)
      gts <- sapply(gts, subset, subset = rep(TRUE, gts.maxcount))
      dimnames(gts) <- list(NULL, gts.rsnumbers)
      
    } else { 
      stop("getGenotypes: Internal error, genotype retrival mode unkown.")
      
    } # end data retrieval, gts is now a character matrix, with snps in columns
    
    if(is.null(gts) || ncol(gts) < 1 || nrow(gts) < 1) {
      warning("getGenotypes: No SNPs left after extraction and formatting of gentoype data")
      return(NULL)    
    }
    
    # replace everything except nucleotides with NA
    gts[grep("[^ACGT]", gts)] <- NA
        
    #############################
    ######## WRITE FILES ########
    #############################

    # this yields a character matrix (converts genotype class columns to its character representation)
    gts <- as.matrix(gts)
    if(is.null(gts) || ncol(gts) < 1 || nrow(gts) < 1)
      return(NULL)
    message(paste("Writing genotypes to ", gwaa.out.fn, " / ", phe.out.fn, " / ", ped.out.fn, " / ", map.out.fn, " ...", sep = ""))
    
    # MAP file
    if(!all(colnames(gts) %in% mapinfo$SNP)) {
      message("Error writing ped/mapfile: Not all genotypes have position information needed to write the mapfile.")
      return(NULL)
    } 
    mapinfo <- mapinfo[match(colnames(gts), mapinfo$SNP), ] # establish order that matches genotype data
    writeLines(
        text = paste(
            mapinfo$CHR, 
            mapinfo$SNP,
            mapinfo$BP, 
            sep = "\t"
        ),
        con = map.out.fn
    )
    # PED file
    # gts can only contain NA or two-character ACGT 
    gts[is.na(gts)] <- "00"
    # when converted from genotype class, are divided by slash
    # gts <- gsub("/", "", gts, fixed = TRUE)
    if(any(nchar(gts) != 2)) {
      warning("getGenotypes: Cannot write genotype file: Illegal genotype matrix format.\n")
      return(NULL)
    }
    write.table(
        x = cbind(
            # header columns, five dummy columns plus one column for individual ID generated by 1:nrow(gts)
            0, 1:nrow(gts), 0, 0, 0, 0,
            # genotypes: introduce a space between two-character genotype representation, as required by ped format
            matrix(data = sub("(.{1})", "\\1 ", gts), ncol = ncol(gts))
        ),
        file = ped.out.fn, 
        sep = "\t", 
        quote = FALSE, 
        row.names = FALSE, 
        col.names = FALSE
    )
    
    # gwaa file (GenABEL)
    capture.output(convert.snp.ped(
        pedfile = ped.out.fn,
        mapfile = map.out.fn, 
        outfile = gwaa.out.fn,
        mapHasHeaderLine = FALSE
    ))
    
    # phe file (GenABEL)
    write.table(
        data.frame(id = 1:nrow(gts), sex = 0, age = 0),
        file = phe.out.fn, 
        row.names = FALSE, 
        sep = "\t"
    )
    
    #######################################################
    ######## TO SNP.DATA OBJECT, REMOVE HOMOZYGOUS ########
    #######################################################
    
    rm(gts)
    capture.output(
      gts <- load.gwaa.data(phenofile = phe.out.fn, genofile = gwaa.out.fn)@gtdata
    )
    
    if(remove.homozygous) {
      maf.snps <- summary(gts)$Q.2
      if(any(maf.snps == 0))      # conditional: avoid subsetting with zero selection (throws error)
        gts <- gts[, maf.snps != 0]
    }
    
  } # end hapmap and linkage retrieval modes, gts is a snp.data object

  if(return.characters) {
    gts <- as.character(gts)
    return(gsub("/", "", gts, fixed = TRUE))
  } else if(return.genotypes) {
    return(as.list(as.genotype(gts)))
  } else {
    return(gts)      
  }

}





parseGtsSourceArg <- function(gts.source) {
  
  res <- list()
  
  # check and parse gts.source argument, set according 'mode' variable for retrieval
  if(!is.null(gts.source) && is.numeric(gts.source) && length(gts.source) == 1) {
    res$mode <- "hapmap"
    res$pop.id <- gts.source
  } else if(!is.null(gts.source) && class(gts.source) == "snp.data" && length(gts.source) == 1) {
    res$mode <- "gwaa-object"
    res$snp.data.obj <- gts.source
  } else if(!is.null(gts.source) && is.character(gts.source) && length(gts.source) == 1 && nchar(gts.source) > 0) {
    if(grepl(".map", gts.source, fixed = TRUE)) {
      res$mode <- "linkage"
      res$map.fn <- gts.source
      res$ped.fn <- sub(".map", ".ped", gts.source, fixed = TRUE)
    } else if(grepl(".ped", gts.source, fixed = TRUE)) {
      res$mode <- "linkage"
      res$map.fn <- sub(".ped", ".map", gts.source, fixed = TRUE)
      res$ped.fn <- gts.source
    } else if(grepl(".gwaa", gts.source, fixed = TRUE)) {
      res$mode <- "gwaa-files"
      res$gwaa.fn <- gts.source
      res$phe.fn <- sub(".gwaa", ".phe", gts.source, fixed = TRUE)
    } else if(grepl(".phe", gts.source, fixed = TRUE)) {
      res$mode <- "gwaa-files"
      res$gwaa.fn <- sub(".phe", ".gwaa", gts.source, fixed = TRUE)
      res$phe.fn <- gts.source
    } else {
      stop("Invalid 'gts.source' argument: Genotype file type unknown, name has to contain .ped, .map, .gwaa or .phe.")
    }
  } else {
    stop("Argument 'gts.source' has to be either a numeric HapMap population ID or a genotype file name or a snp.data object.")    
  }
  
  return(res)
}

