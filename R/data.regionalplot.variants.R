data.regionalplot.variants.single <- function(snps, regions.df, var.options, biomart.config, snpds) {

  id.prune <- var.options$vcf.id.prune
  info.prune <- var.options$vcf.info.prune
  info.color <- var.options$vcf.info.colorize
  af.prune <- var.options$vcf.af.prune # is always set
# das sollte wohl in main...
  if(!is.null(info.prune) && (is.null(names(info.prune)) || length(info.prune) < 1))
    stop("Argument var.options$vcf.info.prune has to be a nonempty, named vector or list.\n")
  if(!is.null(info.color) && (is.null(names(info.color)) || length(info.color) < 1))
    stop("Argument var.options$vcf.info.color has to be a nonempty, named vector or list.\n")
  info.fields <- c(names(info.prune), names(info.color), names(af.prune))
  
  message("Reading vcf file...")
  
  # check that INFO column contains all required fields
  hdr <- scanVcfHeader(var.options$vcf$file)
  if(!all(info.fields %in% rownames(info(hdr))))
    stop(paste("The INFO column of the vcf file does not contain all required fields:", paste(unique(info.fields), collapse = " ; ")))

  # check for matching chromsome names between gwas data and vcf file (scan tabix fails when a chromosome does not exist)
  vcf.chr <- seqnamesTabix(var.options$vcf$file)
  snps.chr <- if(is.null(var.options$vcf$chrom.map)) snps$CHR else merge(snps, var.options$vcf$chrom.map)$CHR.VCF
  for(chr in sort(unique(snps.chr)))
    if(!(chr %in% vcf.chr))
      stop(paste("A query SNP exists on chromosome ", chr, 
              ", but no data exists for that chromosome in the vcf file. Make sure that SNPs for that chromosome exist in the vcf file (region has been resequenced) and /or that chromosome ", 
              chr, " is named accordingly in the vcf file. If not, try to make use of the 'chrom.map' option to remap chromosome names. \nChromosome names found in the vcf file are: ", 
              paste(vcf.chr, collapse = "', '"), ".", sep = "'"
      ))
  
  # all regions as genomic ranges
  range <- GRanges(
    snps.chr, 
    IRanges(start = snps$region.startbp, end = snps$region.endbp, names = as.vector(snps$SNP)), 
    paramRangeID = snps$SNP
  )
  
  if(var.options$vcf$remap.positions == T) {
    # adjust original range by to-remap offset
    range.scan <- range + 2000000
    
    # var.regions = list (regions) of lists (variants in each region)
    # the inner lists can be transformed to data frames with columns according to the VCF file (e.g. CHROM, INFO, GENO.AD.Sample_10893)
    var.regions.scan <- scanVcf(
      var.options$vcf$file, 
      param = ScanVcfParam(fixed = NA, info = NA, geno = NA, which = range.scan)
    )
    names(var.regions.scan) <- as.vector(snps$SNP)
    
    # calculate offset between vcf positions and regions
    offsets <- sapply(
      snps$SNP, 
      function(snp) {
        # var.regions.scan[[snp]]$rowData is a GRanges object
        var.region.df <- data.frame(SNP = names(var.regions.scan[[snp]]$rowData), POS = start(var.regions.scan[[snp]]$rowData))
        match <- merge(var.region.df, regions.df)
        if(nrow(match) < 1) {
          message(paste("Warning: SNPs from region", snp, "were not found as variants in the vcf file (search frame is the region window +- 2 MB - position mapping problem, or SNP identifier problem?). "))
          return(0)
        } else {
          return(sort(match$POS - match$BP.mapped)[nrow(match)/2])
        }
      }
    )
    range <- shift(range, offsets)    
  }

  var.regions <- scanVcf(var.options$vcf$file, param = ScanVcfParam(info = info.fields, which = range))
  
  # make each region a data frame with columns for ID, BP, AF, SPECIAL fields, where SPECIAL is the regex field selected by the user from the INFO data
  var.regions <- lapply(
    var.regions, 
    function(var.region.list) {
      # this code also works for empty var.region.list and returns an empty data frame (with appropriate column names)
      var.region.df <- data.frame(ID = names(var.region.list$rowData))
      var.region.df$BP <- start(var.region.list$rowData)
      for(field in unique(info.fields)) {
#        message(paste("Extracting", field, "field from INFO column..."))
        # in the parsed INFO columns, there may be even multiple elements per field or no data at all, so we have to check before adding to the df
        coldat <- var.region.list$INFO[[field]] # that is a matrix (normally should be single column)
        # column has to have nrow(var.region.df) elements, otherwise it is a nested list
        if(length(unlist(coldat)) == 0) { 
          # need to check for empty coldat because next if-expression fails for that case
          var.region.df[, field] <- unlist(coldat)
        } else if(nrow(var.region.df) == sum(sapply(coldat, function(x) length(unlist(x))))) {
          var.region.df[, field] <- unlist(coldat)
        } else {
          var.region.df[, field] <- unlist(lapply(coldat, paste, collapse = ""))
        }
      }
      # if exist, rename allele frequency field to AF
      if(!is.null(af.prune))
        colnames(var.region.df)[colnames(var.region.df) == names(af.prune)] <- "AF"
      return(var.region.df)
    }
  )
  # the order returnd by scanVcf seems to be currently undetermined. Also, the list names are always by range
  names(range) <- paste(seqnames(range), paste(start(range), end(range), sep = "-"), sep = ":")
  for(idx in 1:length(var.regions)) {
    names(var.regions)[idx] <- elementMetadata(range)[names(range) == names(var.regions)[idx], "paramRangeID"]
  }

  # filter variants and reformat data frame
  var.regions <- lapply(
                  var.regions,
                  function(var.region) {
                    var.region <- vectorElements(var.region)
                    if(nrow(var.region) < 1) return(NULL)
                    # filter
                    if(!is.null(id.prune)) {
                      keep <- rep(FALSE, nrow(var.region))
                      for(regex in id.prune)
                        keep <- keep | grepl(regex, var.region$ID, perl = T)
                      var.region <- var.region[keep, ]
                    }
                    if(!is.null(info.prune)) {
                      keep <- rep(FALSE, nrow(var.region))
                      for(idx in 1:length(info.prune))
                        keep <- keep | grepl(info.prune[idx], var.region[, names(info.prune)[idx]], perl = T)
                    }
                    if(nrow(var.region) < 1) return(NULL)
                    # map to correct position if desired
                    if(var.options$vcf$remap.positions) {
                      # retrieve positions from biomart, DO NOT USE BUFFER HERE (is taken by query SNPs)
                      message("Retrieving variant SNP positions from biomart...")
                      mut.bm <- bm.snps(biomart.config, snpds, as.vector(var.region$ID), use.buffer = FALSE)
                      colnames(mut.bm)[colnames(mut.bm) == "chrom_start"] <- "BP.mapped"
                      # remove multiple chromosomes per SNP as above
                      if(any(duplicated(mut.bm$refsnp_id))) {
                        mut.bm$nchr <- nchar(as.character(mut.bm$chrname))
                        mut.bm <- mut.bm[order(mut.bm$nchr), ]
                        mut.bm <- mut.bm[!duplicated(mut.bm$refsnp_id), ]
                        mut.bm$nchr <- NULL
                      }
                      var.region <- merge(var.region, mut.bm, by.x = "ID", by.y = "refsnp_id", all.x = TRUE)
                      if(all(is.na(var.region$BP.mapped))) {
                        message("Removing variant plot from a region - cannot map positions (try again with remap.positions = FALSE?)")
                        return(NULL)
                      }
                      var.region <- imputeRemappedPositions(var.region)
                    } else {
                      var.region$BP.mapped <- var.region$BP
                    }
                    var.region$COLOR <- rep(hsv(0, 0, 0), nrow(var.region))
                    var.region <- vectorElements(var.region)
                    # switch AF when > 50%
# AF field name?
                    var.region$AF <- var.region$AF * 100
                    var.region[var.region$AF > 50, "AF"] <- 100 - var.region[var.region$AF > 50, "AF"]
                    # set color for variant effect
                    if(!is.null(info.color)) {
                      for(idx in 1:length(info.color)) {
                        var.region.idx <- grep(info.color[idx], var.region[, names(info.color)[idx]], perl = T)
                        if(length(var.region.idx) > 0) 
                          var.region[var.region.idx, "COLOR"] <- rainbow(length(info.color))[idx]
                      }
                    }
                    # data frame has columns ID, BP.mapped, BP, COLOR + fields extracted from info column
                    return(var.region)
                  }
  )
  return(var.regions)
}

data.regionalplot.variants <- function(snps, regions.df, var.options, biomart.config, snpds) {

  if(is.null(var.options$vcf$file)) {

    # dual vcf. Read both files
    var.options1 <- var.options
    var.options2 <- var.options
    var.options1$vcf <- var.options$vcf[[1]]
    var.options2$vcf <- var.options$vcf[[2]]
# TODO assumes that both files use same reference sequence (i.e. reference / minor allele)
# which might not be the case at the moment
    var.regions1 <- data.regionalplot.variants.single(snps, regions.df, var.options1, biomart.config, snpds)
    var.regions2 <- data.regionalplot.variants.single(snps, regions.df, var.options2, biomart.config, snpds)

    var.regions <- sapply(
      names(var.regions1),
      function(idx) {
        if(is.null(var.regions1[[idx]]) & is.null(var.regions2[[idx]])) return(NULL)
        if(is.null(var.regions1[[idx]])) var.regions1[[idx]] <- var.regions2[[idx]][0, ]
        if(is.null(var.regions2[[idx]])) var.regions2[[idx]] <- var.regions1[[idx]][0, ]
        reg.new <- merge(var.regions1[[idx]], var.regions2[[idx]], by = "BP.mapped", all = T)
        reg.new[is.na(reg.new$AF.x), "AF.x"] <- 0
        reg.new[is.na(reg.new$AF.y), "AF.y"] <- 0
        reg.new$COLOR <- reg.new$COLOR.x
        reg.new[is.na(reg.new$COLOR), "COLOR"] <- reg.new[is.na(reg.new$COLOR) ,"COLOR.y"]
        reg.new$ID <- as.character(reg.new$ID.x)
        reg.new[is.na(reg.new$ID), "ID"] <- as.character(reg.new[is.na(reg.new$ID), "ID.y"])
        return(reg.new)
      },
      simplify = F
    )
    # apply af cap
    var.regions <- lapply(var.regions, function(var.region) {
      if(is.null(var.region)) return(NULL)
      var.region <- var.region[var.region$AF.x <= var.options$vcf.af.prune, ]
      var.region <- var.region[var.region$AF.y <= var.options$vcf.af.prune, ]
      if(nrow(var.region) < 1) return(NULL) else var.region
    })

    
  } else {

    var.regions <- data.regionalplot.variants.single(snps, regions.df, var.options, biomart.config, snpds)
    # apply af cap
    var.regions <-lapply(var.regions, function(var.region) {
      if(is.null(var.region)) return(NULL)
      var.region <- var.region[var.region$AF <= var.options$vcf.af.prune, ]
      if(nrow(var.region) < 1) return(NULL) else var.region
    })

  }
  
}
  
