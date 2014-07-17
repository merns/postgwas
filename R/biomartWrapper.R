biomartConfigs <- list(
hsapiens = list(
  snp = list(
    host = "www.biomart.org",  
    mart = "snp", 
    dataset = "hsapiens_snp", 
    filter = list(
      id = "snp_filter"), 
    attr = list(
      id = "refsnp_id", 
      chr = "chr_name", 
      bp = "chrom_start")
  ), 
  gene = list(
    host = "www.biomart.org",
    mart = "ensembl",
    dataset = "hsapiens_gene_ensembl", 
    filter = list(
      id = "entrezgene", 
      name = "hgnc_symbol",
      chr = "chromosome_name", 
      startbp = "start", 
      endbp = "end"), 
    attr = list(
      id = "entrezgene", 
      name = "hgnc_symbol", 
      chr = "chromosome_name", 
      startbp = "start_position", 
      endbp = "end_position", 
      strand = "strand",
      goterms = "name_1006",
      domains = "superfamily")
  ),
  exon = list(
    filter = list(
      geneid = "entrezgene", 
      chr = "chromosome_name"), 
    attr = list(
      chr = "chromosome_name",
      startbp = "exon_chrom_start", 
      endbp = "exon_chrom_end",
      gene.startbp = "start_position",
      gene.endbp = "end_position")
  ), 
  path = list(
      host = "www.biomart.org",
      mart = "REACTOME", 
      dataset = "pathway", 
      filter = list(
          geneid = "referencednasequence_ncbi_id_list"),
# referencednasequence_ensembl_id_list
# referencepeptidesequence_uniprot_id_list
      attr = list(
          geneid = "referencedatabase_ncbi_gene", 
# referencedatabase_ensembl
# referencedatabase_uniprot
          pathid = "stableidentifier_identifier",
          pathname = "_displayname")
  ), 
  proteincomplex = list(
      host = "www.biomart.org",
      mart = "REACTOME", 
      dataset = "complex", 
      filter = list(
          geneid = "referencednasequence_ncbi_id_list"),
      attr = list(
          geneid = "referencedatabase_ncbi_gene",
          complexid = "stableidentifier_identifier",
          complexname = "_displayname") 
  )
), 
athaliana = list(
    snp = list(
        host = "www.biomart.org",
        mart = "plants_variations_16", 
        dataset = "athaliana_eg_snp", 
        filter = list(
            id = "refsnp"), 
        attr = list(
            id = "refsnp_id", 
            chr = "chr_name", 
            bp = "chrom_start")
    ), 
    gene = list(
        host = "www.biomart.org",
        mart = "plants_mart_16",
        dataset = "athaliana_eg_gene", 
        filter = list(
            id = "entrezgene", 
            name = "external_gene_id",
            chr = "chromosome_name", 
            startbp = "start", 
            endbp = "end"), 
        attr = list(
            id = "entrezgene", 
            name = "external_gene_id", 
            chr = "chromosome_name", 
            startbp = "start_position", 
            endbp = "end_position", 
            strand = "strand",
            goterms = "go_name_1006",
            domains = "superfamily_pf")
    ),
    exon = list(
        filter = list(
            geneid = "entrezgene", 
            chr = "chromosome_name"), 
        attr = list(
            chr = "chromosome_name",
            startbp = "exon_chrom_start", 
            endbp = "exon_chrom_end",
            gene.startbp = "start_position",
            gene.endbp = "end_position")
    ), 
    path = list(
        host = "www.biomart.org",
        mart = "REACTOME", 
        dataset = "pathway", 
        filter = list(
            geneid = "referencednasequence_ncbi_id_list"),
        attr = list(
            geneid = "referencedatabase_ncbi_gene", 
            pathid = "stableidentifier_identifier",
            pathname = "_displayname")
    ), 
    proteincomplex = list(
        host = "www.biomart.org",
        mart = "REACTOME", 
        dataset = "complex", 
        filter = list(
            geneid = "referencednasequence_ncbi_id_list"),
        attr = list(
            geneid = "referencedatabase_ncbi_gene",
            complexid = "stableidentifier_identifier",
            complexname = "_displayname") 
    )
),
scerevisiae = list(
  snp = list(
    host = "www.biomart.org",
    mart = "snp", 
    dataset = "scerevisiae_snp", 
    filter = list(
      id = "snp_filter"), 
    attr = list(
      id = "refsnp_id", 
      chr = "chr_name", 
      bp = "chrom_start")
  ), 
  gene = list(
    host = "www.biomart.org",
    mart = "ensembl",
    dataset = "scerevisiae_gene_ensembl", 
    filter = list(
      id = "entrezgene", 
      name = "wikigene_name",
      chr = "chromosome_name", 
      startbp = "start", 
      endbp = "end"), 
    attr = list(
      id = "entrezgene", 
      name = "wikigene_name", 
      chr = "chromosome_name", 
      startbp = "start_position", 
      endbp = "end_position", 
      strand = "strand",
      goterms = "name_1006",
      domains = "superfamily")
  ),
  exon = list(
    filter = list(
      geneid = "entrezgene",
      chr = "chromosome_name"), 
    attr = list(
      chr = "chromosome_name", 
      startbp = "exon_chrom_start", 
      endbp = "exon_chrom_end",
      gene.startbp = "start_position",
      gene.endbp = "end_position")
  ), 
  path = list(
      host = "www.biomart.org",
      mart = "REACTOME", 
      dataset = "pathway", 
      filter = list(
          geneid = "referencednasequence_ncbi_id_list"), 
      attr = list(
          geneid = "referencedatabase_ncbi_gene",
          pathid = "stableidentifier_identifier",
          pathname = "_displayname")
  ), 
  proteincomplex = list(
      host = "www.biomart.org",
      mart = "REACTOME", 
      dataset = "complex", 
      filter = list(
          geneid = "referencednasequence_ncbi_id_list"),
      attr = list(
          geneid = "referencedatabase_ncbi_gene",
          complexid = "stableidentifier_identifier",
          complexname = "_displayname") 
  )
)
)

biomartConfigs$oanatinus <- biomartConfigs$scerevisiae
biomartConfigs$tguttata <- biomartConfigs$scerevisiae
biomartConfigs$fcatus <- biomartConfigs$scerevisiae
biomartConfigs$rnorvegicus <- biomartConfigs$scerevisiae
biomartConfigs$ggallus <- biomartConfigs$scerevisiae
biomartConfigs$ecaballus <- biomartConfigs$scerevisiae
biomartConfigs$pabelii <- biomartConfigs$scerevisiae
biomartConfigs$drerio <- biomartConfigs$scerevisiae
biomartConfigs$tnigroviridis <- biomartConfigs$scerevisiae
biomartConfigs$mdomestica <- biomartConfigs$scerevisiae
biomartConfigs$dmelanogaster <- biomartConfigs$scerevisiae
biomartConfigs$ptroglodytes <- biomartConfigs$scerevisiae
biomartConfigs$sscrofa <- biomartConfigs$scerevisiae
biomartConfigs$mmusculus <- biomartConfigs$scerevisiae
biomartConfigs$btaurus <- biomartConfigs$scerevisiae
biomartConfigs$cfamiliaris <- biomartConfigs$scerevisiae

biomartConfigs$oanatinus$snp$dataset <- "oanatinus_snp"
biomartConfigs$tguttata$snp$dataset <- "tguttata_snp"
biomartConfigs$fcatus$snp$dataset <- "fcatus_snp"
biomartConfigs$rnorvegicus$snp$dataset <- "rnorvegicus_snp"
biomartConfigs$ggallus$snp$dataset <- "ggallus_snp"
biomartConfigs$ecaballus$snp$dataset <- "ecaballus_snp"
biomartConfigs$pabelii$snp$dataset <- "pabelii_snp"
biomartConfigs$drerio$snp$dataset <- "drerio_snp"
biomartConfigs$tnigroviridis$snp$dataset <- "tnigroviridis_snp"
biomartConfigs$mdomestica$snp$dataset <- "mdomestica_snp"
biomartConfigs$dmelanogaster$snp$dataset <- "dmelanogaster_snp"
biomartConfigs$ptroglodytes$snp$dataset <- "ptroglodytes_snp"
biomartConfigs$sscrofa$snp$dataset <- "sscrofa_snp"
biomartConfigs$mmusculus$snp$dataset <- "mmusculus_snp"
biomartConfigs$btaurus$snp$dataset <- "btaurus_snp"
biomartConfigs$cfamiliaris$snp$dataset <- "cfamiliaris_snp"

biomartConfigs$oanatinus$gene$dataset <- "oanatinus_gene_ensembl"
biomartConfigs$tguttata$gene$dataset <- "tguttata_gene_ensembl"
biomartConfigs$hsapiens$gene$dataset <- "hsapiens_gene_ensembl"
biomartConfigs$fcatus$gene$dataset <- "fcatus_gene_ensembl"
biomartConfigs$rnorvegicus$gene$dataset <- "rnorvegicus_gene_ensembl"
biomartConfigs$ggallus$gene$dataset <- "ggallus_gene_ensembl"
biomartConfigs$ecaballus$gene$dataset <- "ecaballus_gene_ensembl"
biomartConfigs$pabelii$gene$dataset <- "pabelii_gene_ensembl"
biomartConfigs$drerio$gene$dataset <- "drerio_gene_ensembl"
biomartConfigs$tnigroviridis$gene$dataset <- "tnigroviridis_gene_ensembl"
biomartConfigs$mdomestica$gene$dataset <- "mdomestica_gene_ensembl"
biomartConfigs$dmelanogaster$gene$dataset <- "dmelanogaster_gene_ensembl"
biomartConfigs$ptroglodytes$gene$dataset <- "ptroglodytes_gene_ensembl"
biomartConfigs$sscrofa$gene$dataset <- "sscrofa_gene_ensembl"
biomartConfigs$mmusculus$gene$dataset <- "mmusculus_gene_ensembl"
biomartConfigs$btaurus$gene$dataset <- "btaurus_gene_ensembl"
biomartConfigs$cfamiliaris$gene$dataset <- "cfamiliaris_gene_ensembl"

biomartConfigs$fcatus$gene$attr$name <- "external_gene_id"
biomartConfigs$rnorvegicus$gene$name <- "rgd_symbol"
biomartConfigs$mmusculus$gene$name <- "mgi_symbol"


# deny retrieval and usage of exons when one of the exon config entries is NA
exons.avail <- function(biomart.config = biomartConfigs$hsapiens, use.buffer = FALSE) {
  # ok when valid buffer exists
  if(use.buffer && !is.null(get("exons.regionalplot", envir = postgwasBuffer))) return(TRUE)
  # deny when attributes are not set
  if(is.null(biomart.config$exon) || any(is.na(unlist(biomart.config$exon)))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# args: a predefined species as character string, e.g. "hsapiens", or a list with all required biomart field names. Required fields are listed by calling names(biomartConfigs$hsapiens)
# value: a biomaRt dataset object
bm.init.snps <- function(species = NULL) {

  if(is.null(species))
    stop(paste("The following predefined configs are available (use as quoted string in the function call):", paste(names(biomartConfigs), collapse = ", ")))

  if(is.list(species)) { # user-defined biomart attribute list
    config <- species
  } else if(is.null(biomartConfigs[[species]])) {
      stop(paste("Biomart config not available for selected species. Available predefined configs:", paste(names(biomartConfigs), collapse = ", ")))
  } else {
    config <- biomartConfigs[[species]]
  }

  message("Initializing biomart connection (SNPs)...")
  ds <- useMart(biomart = config$snp$mart, dataset = config$snp$dataset, host = config$snp$host)

  # check attribs exist
  if(sum(listAttributes(ds)[, 1] %in% unlist(config$snp$attr)) < length(config$snp$attr)) 
    stop(paste("One of the biomart attributes [", paste(config$snp$attr, collapse = "; "), "] in dataset [", config$snp$dataset, "] does not exist"))
  # filters
  if(sum(listFilters(ds)[, 1] %in% unlist(config$snp$filter)) < length(config$snp$filter)) 
    stop(paste("One of the biomart filters [", paste(config$snp$filter, collapse = "; "), "] in dataset [", config$snp$dataset, "] does not exist"))

  return(ds)
}

# args: a predefined species as character string, e.g. "hsapiens", or a list with all required biomart field names. Required fields are listed by calling names(biomartConfigs$hsapiens)
# value: a biomaRt dataset object
bm.init.genes <- function(species = NULL) {

  if(is.null(species))
    stop(paste("The following predefined configs are available (use as quoted string in the function call):", paste(names(biomartConfigs), collapse = ", ")))

  if(is.list(species)) { # user-defined biomart attribute list
    config <- species
  } else if(is.null(biomartConfigs[[species]])) {
      stop(paste("Biomart config not available for selected species. Available predefined configs:", paste(names(biomartConfigs), collapse = ", ")))
  } else {
    config <- biomartConfigs[[species]]
  }

  message("Initializing biomart connection (genes)...")
  ds <- useMart(biomart = config$gene$mart, dataset = config$gene$dataset, host = config$gene$host)

  # check attribs exist
  if(sum(listAttributes(ds)[, 1] %in% unlist(config$gene$attr)) < length(config$gene$attr)) 
    stop(paste("One of the biomart attributes", paste(config$gene$attr, collapse = "; "), "in dataset", config$gene$dataset, "does not exist"))
  # exon attribs
  if(exons.avail(config, FALSE) && sum(listAttributes(ds)[, 1] %in% unlist(config$exon$attr)) < length(config$exon$attr)) 
    stop(paste("One of the biomart attributes", paste(config$exon$attr, collapse = "; "), "in dataset", config$gene$dataset, "does not exist"))
  # gene filters
  if(sum(listFilters(ds)[, 1] %in% unlist(config$gene$filter)) < length(config$gene$filter)) 
    stop(paste("One of the biomart filters", paste(config$gene$filter, collapse = "; "), "in dataset", config$gene$dataset, "does not exist"))

  return(ds)
}


# args: a predefined species as character string, e.g. "hsapiens", or a list with all required biomart field names. Required fields are listed by calling names(biomartConfigs$hsapiens)
# value: a biomaRt dataset object
bm.init.path <- function(species = NULL) {
  
  if(is.null(species))
    stop(paste("The following predefined configs are available (use as quoted string in the function call):", paste(names(biomartConfigs), collapse = ", ")))
  
  if(is.list(species)) { # user-defined biomart attribute list
    config <- species
  } else if(is.null(biomartConfigs[[species]])) {
    stop(paste("Biomart config not available for selected species. Available predefined configs:", paste(names(biomartConfigs), collapse = ", ")))
  } else {
    config <- biomartConfigs[[species]]
  }
  
  message("Initializing biomart connection (pathways)...")
  ds <- useMart(biomart = config$path$mart, dataset = config$path$dataset, host = config$path$host)
  
  # check attribs exist
  if(sum(listAttributes(ds)[, 1] %in% unlist(config$path$attr)) < length(config$path$attr)) 
    stop(paste("One of the biomart attributes", paste(config$path$attr, collapse = "; "), "in dataset", config$path$dataset, "does not exist"))
  # filters
  if(sum(listFilters(ds)[, 1] %in% unlist(config$path$filter)) < length(config$path$filter)) 
    stop(paste("One of the biomart filters", paste(config$path$filter, collapse = "; "), "in dataset", config$path$dataset, "does not exist"))
  
  return(ds)
}



# args: a predefined species as character string, e.g. "hsapiens", or a list with all required biomart field names. Required fields are listed by calling names(biomartConfigs$hsapiens)
# value: a biomaRt dataset object
bm.init.proteincomplex <- function(species = NULL) {
  
  if(is.null(species))
    stop(paste("The following predefined configs are available (use as quoted string in the function call):", paste(names(biomartConfigs), collapse = ", ")))
  
  if(is.list(species)) { # user-defined biomart attribute list
    config <- species
  } else if(is.null(biomartConfigs[[species]])) {
    stop(paste("Biomart config not available for selected species. Available predefined configs:", paste(names(biomartConfigs), collapse = ", ")))
  } else {
    config <- biomartConfigs[[species]]
  }
  
  message("Initializing biomart connection (protein interaction)...")
  ds <- useMart(biomart = config$proteincomplex$mart, dataset = config$proteincomplex$dataset, host = config$proteincomplex$host)
  
  # check attribs exist
  if(sum(listAttributes(ds)[, 1] %in% unlist(config$proteincomplex$attr)) < length(config$proteincomplex$attr)) 
    stop(paste("One of the biomart attributes", paste(config$proteincomplex$attr, collapse = "; "), "in dataset", config$proteincomplex$dataset, "does not exist"))
  # filters
  if(sum(listFilters(ds)[, 1] %in% unlist(config$proteincomplex$filter)) < length(config$proteincomplex$filter)) 
    stop(paste("One of the biomart filters", paste(config$proteincomplex$filter, collapse = "; "), "in dataset", config$proteincomplex$dataset, "does not exist"))
  
  return(ds)
}


# returns a data frame with the columns 
# "refsnp_id", "chrname", "chrom_start" 
# in the listed order. Columns are cast to vector. 
bm.snps <- function(
              config = biomartConfigs$hsapiens, 
              ds = bm.init.snps(config), 
              filter.val, 
              use.buffer = FALSE,
              remove.dupl = TRUE) {
  attrnames <- c(config$snp$attr$id, config$snp$attr$chr, config$snp$attr$bp)
  if(use.buffer && !is.null(get("snps", envir = postgwasBuffer))) {
    message("Using buffer data (SNPs)")
    buf <- get("snps", envir = postgwasBuffer)
    if(!all(attrnames == colnames(buf)))
      stop("Error in snp buffer data - clear buffer (run clearPostgwasBuffer()) or set use.buffer = FALSE")
    else
      snps <- buf[buf[, config$snp$attr$id] %in% filter.val, ]
  } else {
    snps <- getBM(
              attributes = attrnames, 
              filters = config$snp$filter$id, 
              values = filter.val, 
              mart = ds
            )
    if(use.buffer) {
      assign("snps", snps, envir = postgwasBuffer)
    }
  }
  colnames(snps) <- c("refsnp_id", "chrname", "chrom_start")
  
  # rare cases of multiple positions (chromosomes) per SNP in biomart can be removed (e.g. MHC assemblies)
  # keep the smallest chromsome name by character length
  # duplicate removes elements with larger index -> sort by character length
  if(remove.dupl) 
    snps <- snps[!(snps$refsnp_id %in% snps$refsnp_id[duplicated(snps$refsnp_id)]), ]

  snps$chrname <- as.vector(snps$chrname)
  return(snps)
}


# returns a data frame with the columns 
# "geneid", genename", "chrname", "start", "end"
# in the listed order.
# We can have strange chromosomes for genes like LRG genes, but this should not matter
bm.genes <- function(
                        config = biomartConfigs$hsapiens, 
                        ds = bm.init.genes(config), 
                        use.buffer = FALSE) {
  attrnames <- c(config$gene$attr$id, config$gene$attr$name, config$gene$attr$chr, config$gene$attr$startbp, config$gene$attr$endbp)
  if(use.buffer && !is.null(get("genes", envir = postgwasBuffer))) {
    message("Using buffer data (genes)")
    genes <- get("genes", envir = postgwasBuffer)
    if(!all(attrnames %in% colnames(genes)))
      stop("Error in genes buffer data - clear buffer (run clearPostgwasBuffer()) or set use.buffer = FALSE")
  } else {
    genes <- getBM(attributes = attrnames, mart = ds)
    if(use.buffer) {
      assign("genes", genes, envir = postgwasBuffer)
    }
  }
  names(genes) <- c("geneid", "genename", "chrname", "start", "end")
  genes$geneid <- as.vector(genes$geneid)
  genes$genename <- as.vector(genes$genename)
  genes$chrname <- as.vector(genes$chrname)
  return(genes)
}


# returns a data frame with the columns 
# id, name, chr, startbp, endbp, strand 
# in the listed order.
bm.genes.regionalplot <- function(
                            config = biomartConfigs$hsapiens, 
                            ds = bm.init.genes(config), 
                            chr,      # position can be vectors for multiple regions
                            bp.start, 
                            bp.end, 
                            use.buffer = FALSE) {
  attrnames <- c(config$gene$attr$id, config$gene$attr$name, config$gene$attr$chr, config$gene$attr$startbp, config$gene$attr$endbp, config$gene$attr$strand)
  if(use.buffer && !is.null(get("genes.regionalplot", envir = postgwasBuffer))) {
    message("Using buffer data (genes)")
    genes <- get("genes.regionalplot", envir = postgwasBuffer)
    if(!all(attrnames %in% colnames(genes)))
      stop("Error in genes buffer data - clear buffer (run clearPostgwasBuffer()) or set use.buffer = FALSE")
  } else {
    genes <- unique(list2df(lapply(
      1:length(chr), 
      function(idx) 
        getBM( 
          attributes = attrnames, 
          filters = c(config$gene$filter$chr, config$gene$filter$startbp, config$gene$filter$endbp),
          values = list(chr[idx], bp.start[idx], bp.end[idx]),
          mart = ds 
        )
    )))
    if(use.buffer) {
      assign("genes.regionalplot", genes, envir = postgwasBuffer)
    }
  }
  names(genes) <- c("id", "name", "chromo", "start", "end", "strand")
  genes$id <- as.vector(genes$id)
  genes$name <- as.vector(genes$name)
  genes$mean <- genes$start + abs((genes$start - genes$end)/2)
  return(genes)
}


# using the buffer, this function returns all exons in the buffer!
# filtering for chromosome is mandatory because highly variable regions like MHC can have multiple chromos for the same position
bm.exons <- function(
                   config = biomartConfigs$hsapiens, 
                   ds = bm.init.genes(config), 
                   filter.val = list(geneid = c(""), chromo = c("")), 
                   use.buffer = FALSE) {
  names(filter.val) <- c(config$exon$filter$geneid, config$exon$filter$chr)
  attrnames <- c(config$exon$attr$chr, config$exon$attr$startbp, config$exon$attr$endbp, config$exon$attr$gene.startbp, config$exon$attr$gene.endbp)
  if(use.buffer && !is.null(get("exons.regionalplot", envir = postgwasBuffer))) {
    message("Using buffer data (exons)")
    exons <- get("exons.regionalplot", envir = postgwasBuffer)
    if(!all(attrnames == colnames(exons)))
      stop("Error in exon buffer data - clear buffer (run clearPostgwasBuffer()) or set use.buffer = FALSE")
  } else {
    exons <- getBM(
      attributes = attrnames, 
      filters = c(config$exon$filter$geneid, config$exon$filter$chr), 
      values = filter.val, 
      mart = ds
    )
    if(use.buffer) {
      assign("exons.regionalplot", exons, envir = postgwasBuffer)
    }
  }
  names(exons) <- c("chromo", "start", "end", "genestart", "geneend")
  return(exons)
}

bm.genes.goterms <- function(
                       config = biomartConfigs$hsapiens, 
                       ds = bm.init.genes(config),
                       filter.name = biomartConfigs$hsapiens$gene$filter$name, 
                       filter.val, 
                       use.buffer = FALSE
                     ) {
  
  if(is.null(config$gene$attr$goterms)) 
    stop("Biomart configuration is incomplete - the 'goterms' element has to be specified for the current operation\n")
  
  attrnames <- c(config$gene$attr$id, config$gene$attr$name, config$gene$attr$goterms)
  
  if(use.buffer && !is.null(get("goterms", envir = postgwasBuffer))) {
    message("Using buffer data (GO terms)")
    goterms <- get("goterms", envir = postgwasBuffer)
    if(!all(attrnames == colnames(goterms)))
      stop("Error in GO term buffer data - clear buffer (run clearPostgwasBuffer()) or set use.buffer = FALSE")
  } else {
    goterms <- getBM( 
        attributes = attrnames, 
        filters = filter.name,
        values = filter.val,
        mart = ds 
    )
    if(use.buffer) {
      assign("goterms", goterms, envir = postgwasBuffer)
    }
  }
  
  colnames(goterms) <- c("geneid", "genename", "goterm.name")
  return(goterms)
}

bm.genes.domains <- function(
                       config = biomartConfigs$hsapiens, 
                       ds = bm.init.genes(config),
                       filter.name = biomartConfigs$hsapiens$gene$filter$name, 
                       filter.val
                     ) {
  
  if(is.null(config$gene$attr$domains)) 
    stop("Biomart configuration is incomplete - the 'domains' element has to be specified for the current operation\n")
  
  attrnames <- c(config$gene$attr$id, config$gene$attr$name, config$gene$attr$domains)
  domains <- getBM( 
               attributes = attrnames, 
               filters = filter.name,
               values = filter.val,
               mart = ds 
             )
  colnames(domains) <- c("geneid", "genename", "domain.name")
  domains <- domains[domains$domain.name != "", ]
  return(domains)
}


bm.id2name <- function(
                filter.ids = NULL,
                config = biomartConfigs$hsapiens, 
                ds = bm.init.genes(config), 
                use.buffer = FALSE
              ) {

  if(is.null(filter.ids))
    stop("Argument filter.ids has to be specified.\n")
  
  if(use.buffer && !is.null(get("genes", envir = postgwasBuffer))) {
    message("Using buffer data (genes)")
    buf <- get("genes", envir = postgwasBuffer)
    if(!all(c(config$gene$attr$id, config$gene$attr$name) %in% colnames(buf)))
      stop("Error in genes buffer data - clear buffer (run clearPostgwasBuffer()) or set use.buffer = FALSE")
  } else {
    use.buffer <- FALSE
  }
  
  if(is.data.frame(filter.ids) & config$gene$attr$id %in% colnames(filter.ids)) {
    if(use.buffer) {
      map <- merge(filter.ids, buf, all.x = TRUE)
      map <- map[, c(colnames(filter.ids), config$gene$attr$name)]
    } else {
      map <- getBM(
        attributes = c(config$gene$attr$id, config$gene$attr$name), 
        filters = config$gene$filter$id,
        values = filter.ids[, config$gene$attr$id],
        mart = ds 
      )
      map <- merge(filter.ids, map, all.x = TRUE)
    }
  } else {
    if(use.buffer) {
      map <- buf[
               buf[, config$gene$attr$id] %in% as.vector(filter.ids), 
               c(config$gene$attr$id, config$gene$attr$name)
             ]
    } else {
      map <- getBM(
        attributes = c(config$gene$attr$id, config$gene$attr$name), 
        filters = config$gene$filter$id,
        values = filter.ids,
        mart = ds 
      )
    }
  }
  # make names unique, keep NAs
  map <- map[!duplicated(map[, config$gene$attr$name]) | is.na(map[, config$gene$attr$name]), ]
  return(map)

}


bm.name2id <- function(
                filter.names = NULL,
                config = biomartConfigs$hsapiens, 
                ds = bm.init.genes(config), 
                use.buffer = FALSE
              ) {
  if(is.null(filter.names))
    stop("Argument filter.names has to be specified.\n")
  
  if(use.buffer && !is.null(get("genes", envir = postgwasBuffer))) {
    message("Using buffer data (genes)")
    buf <- get("genes", envir = postgwasBuffer)
    if(!all(c(config$gene$attr$id, config$gene$attr$name) %in% colnames(buf)))
      stop("Error in genes buffer data - clear buffer (run clearPostgwasBuffer()) or set use.buffer = FALSE")
  } else {
    use.buffer <- FALSE
  }
  
  if(is.data.frame(filter.names) & config$gene$attr$name %in% colnames(filter.names)) {
    if(use.buffer) {
      map <- merge(filter.names, buf, all.x = TRUE)
      map <- map[, c(colnames(filter.names), config$gene$attr$id)]
    } else {
      map <- getBM(
        attributes = c(config$gene$attr$id, config$gene$attr$name), 
        filters = config$gene$filter$name,
        values = filter.names[, config$gene$attr$name],
        mart = ds 
        )
      map <- merge(filter.names, map, all.x = TRUE)
    }
  } else {
    if(use.buffer) {
      map <- buf[
          buf[, config$gene$attr$name] %in% as.vector(filter.names), 
          c(config$gene$attr$id, config$gene$attr$name)
        ]
    } else {
      map <- getBM(
        attributes = c(config$gene$attr$id, config$gene$attr$name), 
        filters = config$gene$filter$name,
        values = filter.names,
        mart = ds 
      )
    }
  }
  # make IDs unique, keep NAs
  map <- map[!duplicated(map[, config$gene$attr$id]) | is.na(map[, config$gene$attr$id]), ]
  return(map)
}


bm.remapSnps <- function(
                  file, 
                  snp.col = "SNP", 
                  config = biomartConfigs$hsapiens, 
                  ds = bm.init.snps(config), 
                  use.buffer = FALSE, 
                  toFile = paste(file, ".remapped", sep = "")
                ) {
  
  if(is.numeric(snp.col)) {
    f <- na.omit(read.table(file, header=FALSE))
    colnames(f)[snp.col] <- "SNP" 
    snp.col <- "SNP"
  } else {
    f <- na.omit(read.table(file, header=TRUE))
  }
  
  if(!(snp.col %in% colnames(f)))
    stop(paste("Source file does not contain a column '", snp.col, "' or is not readable by the read.table function.\n"))
  
  message("Retrieving SNP positions from biomart, this can take a while.")
  snps.bm <- bm.snps(
    filter.val = f[, snp.col], 
    config = config, 
    use.buffer = use.buffer
  )
  
  colnames(snps.bm)[colnames(snps.bm) %in% c("chrname", "chrom_start")] <- c("CHR", "BP")
  f.new <- merge(f, snps.bm, by.x = snp.col, by.y = "refsnp_id", suffixes = c(".original", ""), all.x = TRUE)
  f.new <- f.new[match(f[, snp.col], f.new[, snp.col]), ] # restore original order
  
  if(!is.null(toFile) && is.character(toFile) && length(toFile) == 1)
    write.table(f.new, toFile, row.names = F, col.names = T)
  
}
