## returns a data frame "geneid.x", "geneid.y", "label"
#getInteractions <- function(
#    filter.ids,
#    biomart.config = biomartConfigs$hsapiens,
#    path = TRUE, 
#    ppi = TRUE, 
#    domains = TRUE, 
#    GO = FALSE,
#    GOpackagename = "org.Hs.eg.db",
#    toFile = "postgwas.interaction.download"
#) {
#
#  if(missing(filter.ids) || is.null(filter.ids) || !is.vector(filter.ids))
#    stop("Argument 'filter.ids' has to be a non-NULL vector.\n")
#  filter.ids <- as.vector(filter.ids)
#  
#  net.all <- data.frame(geneid.x = 1, geneid.y = 1, label = 1)[-1, ] # empty data frame
#  
#  if(path) {
#    net <- getInteractions.path(filter.ids = filter.ids, biomart.config = biomart.config, toFile = NULL)
#    net$label <- "path"
#    net.all <- rbind(net.all, net)
#  }
#  
#  if(ppi) {
#    net <- getInteractions.proteincomplex(filter.ids = filter.ids, biomart.config = biomart.config, toFile = NULL)
#    net$label <- "ppi"
#    net.all <- rbind(net.all, net)
#  }
#    
#  if(domains) {
#    net <- getInteractions.domains(filter.ids = filter.ids, biomart.config = biomart.config, filter.type = biomart.config$gene$filter$id, toFile = NULL)
#    net$label <- "domain"
#    net.all <- rbind(net.all, net[, c("geneid.x", "geneid.y", "label")])
#  }
#  
#  if(GO) {
#    message("Detecting shared GO term architecture...")
#    if(GOpackagename)
#      # returns a data frame "geneid.x", "geneid.y", "weight"
#      net <- getInteractions.GO(
#          filter.ids = filter.ids, 
#          toFile = NULL
#      )
#    net$label <- "GO"
#    net.all <- rbind(net.all, net[, c("geneid.x", "geneid.y", "label")])
#  }
#  
#  if(!is.null(toFile) && toFile != "")
#    write.table(net.all, toFile, sep = "\t", col.names = TRUE, row.names = FALSE)
#  
#  return(net.all)
#}



# returns a data frame "geneid.x", "geneid.y", "label"
getInteractions.path <- function(
    filter.ids, 
    biomart.config = biomartConfigs$hsapiens,
    ds = bm.init.path(biomart.config), 
    toFile = "postgwas.interaction.download"
) {
  
  if(missing(filter.ids) || is.null(filter.ids) || !is.vector(filter.ids))
    stop("Argument 'filter.ids' has to be a non-NULL vector.\n")
  filter.ids <- as.vector(filter.ids)
  
  dat <- getBM(
      attributes = c(biomart.config$path$attr$pathid, biomart.config$path$attr$pathname, biomart.config$path$attr$geneid), 
      filters = biomart.config$path$filter$geneid, 
      values = filter.ids, 
      mart = ds
  )
  colnames(dat)[3] <- c("geneid")
  
  message("Constructing network...")
  mdat <- merge(dat, dat, by = biomart.config$path$attr$pathid)
  
  # proper column order and labels for gwas2network
  network <- mdat[, c("geneid.x", "geneid.y", paste(biomart.config$path$attr$pathname, "x", sep = "."))]
  colnames(network)[3] <- "label"
  
  if(!is.null(toFile) && toFile != "")
    write.table(network, toFile, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  return(network)
}



# returns a data frame "geneid.x", "geneid.y", "label"
getInteractions.proteincomplex <- function(
    filter.ids, 
    biomart.config = biomartConfigs$hsapiens,
    ds = bm.init.proteincomplex(biomart.config), 
    toFile = "postgwas.interaction.download"
) {
  
  if(missing(filter.ids) || is.null(filter.ids) || !is.vector(filter.ids))
    stop("Argument 'filter.ids' has to be a non-NULL vector.\n")
  filter.ids <- as.vector(filter.ids)
  
  dat <- getBM(
      attributes = c(biomart.config$proteincomplex$attr$complexid, biomart.config$proteincomplex$attr$complexname, biomart.config$proteincomplex$attr$geneid), 
      filters = biomart.config$proteincomplex$filter$geneid, 
      values = filter.ids, 
      mart = ds
  )
  colnames(dat)[3] <- c("geneid")
  
  message("Constructing network...")
  mdat <- merge(dat, dat, by = biomart.config$proteincomplex$attr$complexid)
  
  # proper column order and labels for gwas2network
  network <- mdat[, c("geneid.x", "geneid.y", paste(biomart.config$proteincomplex$attr$complexname, "x", sep = "."))]
  colnames(network)[3] <- "label"
  
  if(!is.null(toFile) && toFile != "")
    write.table(network, toFile, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  return(network)
}

# legacy biogrid download code
#getInteractions.ppi <- function(
#                         filter.ids = c(), 
#                         additionalIdentifierTypes = c("OFFICIAL_SYMBOL", "ENTREZ_GENE", "ENSEMBL", "SWISSPROT"), 
#                         taxId = 9606,
#                         includeInteractors = FALSE, 
#                         sourceDatabaseList = c("BioGRID", "INTACT", "MINT", "DIP"), 
#                         toFile = "postgwas.interaction.download"
#                   ) {
#
#  if(!is.null(filter.ids) && is.vector(filter.ids) && length(filter.ids) > 0) {
#    filter.ids <- unique(as.vector(filter.ids))
#    if(is.null(additionalIdentifierTypes))
#      stop("Argument 'additionalIdentifierTypes' has to be set\n")
#    if(is.null(includeInteractors))
#      stop("Argument 'includeInteractors' has to be set\n")
#    filterquery <- paste(
#                     "geneList=", 
#                     paste(filter.ids, collapse = "|"), "&", 
#                     "additionalIdentifierTypes=",
#                     paste(additionalIdentifierTypes, collapse = "|"), "&",
#                     "includeInteractorInteractions=FALSE", "&", 
#                     "includeInteractors=", includeInteractors,
#                     sep = ""
#                   )
#    # for many filter ids, we have to divide the queries (URL length cap)
#    # when we do that, including interactor interactions is invalid (this feature is disabled anyways)
#    if(length(filter.ids) > 500) {
#      res <- list()
#      filter.ids.part <- lapply(
#                           seq(1, length(filter.ids), by = 500),  # the split points
#                           function(split)
#                             as.vector(na.omit(filter.ids[split:(split+499)]))
#                         )
#      res <- lapply(
#               filter.ids.part, 
#               function(part) 
#                 getInteractions.ppi(
#                   filter.ids = part,
#                   taxId = taxId,
#                   additionalIdentifierTypes = additionalIdentifierTypes, 
#                   includeInteractors = includeInteractors, 
#                   sourceDatabaseList = sourceDatabaseList, 
#                   toFile = NULL
#                 ) 
#             )
#      return(unique(list2df(res)))
#    }
#  } else {
#    warning("Argument 'filter.ids' is not a vector with length > 0. Downloading all interactions.\n")
#    filterquery = NULL
#  }
#
#  # determine number of results to retrieve - there is a cap of 10000 per query
#  message("Determining the number of interactions to download...")
#  url.rows <- url(paste(
#                  "http://webservice.thebiogrid.org/resources/interactions?", 
#                  "selfInteractionsExcluded=true&format=count&", 
#                  filterquery, "&", 
#                  "taxId=",
#                  taxId,
#                  sep = ""
#                ))
#
#  rows <- as.numeric(readLines(url.rows, warn = FALSE))
#  close(url.rows)
#  
#  if(rows == 0)
#    return(data.frame(genename.x = "", genename.y = "")[-1, ])
#  
#  queries.data <- paste(
#                    "http://webservice.thebiogrid.org/resources/interactions?",
#                    "selfInteractionsExcluded=true&format=tab1&",
#                    filterquery, "&",
#                    "sourceDatabaseList=",
#                    paste(sourceDatabaseList, collapse = "|"), 
#                    "&taxId=",
#                    taxId, 
#                    "&start=",
#                    gsub(" ", "", format(0:(rows/10000) * 10000, scientific = FALSE)), 
#                    sep = ""
#                  )
#  queries.data[1] <- paste(queries.data[1], "includeHeader=true", sep = "&")
#  urls.data <- lapply(queries.data, url)
#
#  # paste together the repeated retrieves and download
#  message(paste("Downloading", rows, "interactions"))
#  res <- lapply(
#           urls.data, 
#           function(url) {
#             message(".")
#             dl <- list2df(strsplit(readLines(url, warn = FALSE), split = "\t", fixed = TRUE))
#             close(url)
#             return(dl)
#           }
#         )
#  message("")
#  res <- as.data.frame(list2df(res))
#  colnames(res) <- as.matrix(res)[1, ]
#  res <- res[-1, ]  #remove header
#
#  if(!is.null(toFile) && toFile != "")
#    write.table(res, toFile, sep = "\t", col.names = TRUE, row.names = FALSE)
#
#  return(
#    data.frame(
#      genename.x = res[, "OFFICIAL_SYMBOL_A"], 
#      genename.y = res[, "OFFICIAL_SYMBOL_B"], 
#      stringsAsFactors = FALSE
#    )
#  )
#}



# returns a data frame "genename.x", "genename.y", "geneid.x", "geneid.y", "domain.name"
getInteractions.domains <- function(
    filter.ids, 
    biomart.config = biomartConfigs$hsapiens,
    filter.type = biomartConfigs$hsapiens$gene$filter$name,
    ds = bm.init.genes(biomart.config),
    min.occurence = NULL, 
    max.occurence = if(length(filter.ids) < 20) 
                      NULL 
                    else 
                      length(filter.ids) / (logb(length(filter.ids), 3) +1),
    toFile = "postgwas.interaction.download"
) {
  
  if(missing(filter.ids) || is.null(filter.ids) || !is.vector(filter.ids))
    stop("Argument 'filter.ids' has to be a non-NULL vector.\n")
  
  filter.ids <- as.vector(filter.ids)
  
  if(is.null(filter.type))
    stop("Argument 'filter.type' has to be set.\n")
  
  dom <- bm.genes.domains(
           config = biomart.config, 
           ds = ds, 
           filter.name = filter.type, 
           filter.val = filter.ids
         )

  dom.net <- merge(dom, dom, by = "domain.name")
  
  if(!is.null(min.occurence)) {
    keep <- table(dom.net$domain.name)
    keep <- names(keep)[keep > min.occurence]
    dom.net <- dom.net[dom.net$domain.name %in% keep, ]  
  }
  
  if(!is.null(max.occurence)) {
    keep <- table(dom.net$domain.name)
    keep <- names(keep)[keep < max.occurence]
    dom.net <- dom.net[dom.net$domain.name %in% keep, ]  
  }
  
  dom.net <- dom.net[, c(3,5,2,4,1)]
  
  if(!is.null(toFile) && toFile != "")
    write.table(dom.net, toFile, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  return(dom.net)
}




# returns a data frame "geneid.x", "geneid.y", "weight"
utils::globalVariables("GOSimEnv")
getInteractions.GO <- function(
                             filter.ids, 
                             GOpackagename = "org.Hs.eg.db", 
                             ontology = "BP", 
                             toFile = "postgwas.interaction.download", 
                             cores = 1, 
                             ...
                           ) {

  # these packages are required for GOSim functionality (see GoSimRedist.R)
  for(pkg in c("GO.db")) {
    if(!(pkg %in% names(installed.packages()[, 'Package']))) {
      if(interactive() && readline(paste("Package '", pkg, "' is not installed. Try to install? [Y/N] ")) %in% c("Y", "y")) {
        biocLite <- NULL
        source("http://bioconductor.org/biocLite.R")
        biocLite(pkg)
      } else {
        stop(paste("This function requires the package '", pkg, "' to be installed.\n"))		
      }
    }
    suppressPackageStartupMessages(do.call(require, list(package = pkg, quietly = TRUE)))
  }

  if(missing(filter.ids) || is.null(filter.ids) || !is.vector(filter.ids) || length(filter.ids) <= 0)
    stop("Argument 'filter.ids' has to be a non-NULL vector of entrez geneids.\n")
  
  # load proper GO organism data
  success <- do.call(library, list(package = GOpackagename, logical.return = TRUE))
  if(!success)
	  stop("The GeneOnotology annotation package supplied in argument 'GOpackagename' is not installed or cannot be loaded.\n")
  
  # init GOSim with the selected GO organism annotation package (and given ontology)
  org.name.new <- get(grep("ORGANISM", ls(paste("package", GOpackagename, sep = ":")), value = TRUE)[1])
  org.dat.new <- get(grep("GO", ls(paste("package", GOpackagename, sep = ":")), value = TRUE)[1])
  if(is.null(GOSimEnv$organism) || sub("human", "Homo sapiens", GOSimEnv$organism) != org.name.new) {
    setEvidenceLevel(organism = org.name.new, gomap = org.dat.new, cores = cores)
    setOntology(ontology, loadIC = FALSE, cores = cores)
  } else {
    if(is.null(GOSimEnv$ontology) || GOSimEnv$ontology != ontology) 
      setOntology(ontology, loadIC = FALSE, cores = cores)
  }
  
  message("Starting gene similarity calculation with GOSim (this can take a while)...")
  simmatrix <- getGeneSim(as.vector(filter.ids),cores = cores, ...)
  
  res <- list2df(lapply(
    colnames(simmatrix), 
    function(colname) 
      data.frame(
        geneid.x = colname, 
        geneid.y = rownames(simmatrix), 
        weight = simmatrix[, colname]
        )
    ))

  res <- res[res$weight > 0, ]
  
  if(!is.null(toFile) && toFile != "")
    write.table(res, toFile, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  return(res)
  
}

