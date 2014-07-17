gwasGOenrich <- function(
    gwas, 
    ontology = "BP", 
    pkgname.GO = "org.Hs.eg.db",
    topGOalgorithm = "weight01",
    pruneTermsBySize = 10,
    plotSigTermsToFile = 5
) {
  
  # check gwas argument
  if(
      missing(gwas) || 
      is.null(gwas) || 
      !is.data.frame(gwas) || 
      !("gene.p" %in% colnames(gwas)) || 
      !any(c("geneid", "genename") %in% colnames(gwas)) || 
      nrow(gwas) < 1
  )
    stop("Argument 'gwas' has to be a data frame containing the columns 'gene.p' and one of 'geneid' or 'genename' and at least one row of data.")
  
  # check pruneTermsBySize argument
  if(is.null(pruneTermsBySize) || !is.numeric(pruneTermsBySize) || pruneTermsBySize < 1) {
    pruneTermsBySize <- 1
    warning("Invalid argument 'pruneTermsBySize' has been set to 1 (no pruning).")
  }
  
  # check ontology arg
  if(is.null(ontology) || !is.character(ontology) || length(ontology) > 1 || !(ontology %in% c("BP", "MF", "CC")))
    stop("Argument 'ontology' can only be 'BP', 'MF' or 'CC'.")

  # check pkgname.GO argument and load pkgname.GO library
  if(is.null(pkgname.GO) || !is.character(pkgname.GO) || !(pkgname.GO %in% installed.packages()))
    stop("Argument 'pkgname.GO' has to be a character string matching the name of an installed GO annotation package.")
  capture.output(suppressPackageStartupMessages(require(pkgname.GO, character.only = TRUE)))
  
  suppressPackageStartupMessages(require(topGO))
  
  # check topGOalgorithm arg
  if(!(topGOalgorithm %in% whichAlgorithms())) 
    stop(paste("Argument 'topGOalgorithm' can only be one of: ", paste(whichAlgorithms(), collapse = ", ")))
  
  # auto-detect identifier used
  geneid.type <- NULL
  message("Auto-detect gene identifier type used... ", appendLF = FALSE)
  if("geneid" %in% colnames(gwas)) {
    if(sum(unique(as.vector(gwas$geneid)) %in%  names(unlist(AnnotationDbi::as.list(get(paste(substr(pkgname.GO, 1, nchar(pkgname.GO) -3), "SYMBOL", sep = "")))))) > length(unique(gwas$geneid)) / 2)
      geneid.type <- "entrez"
    else if(sum(unique(as.vector(gwas$geneid)) %in%  unlist(AnnotationDbi::as.list(get(paste(substr(pkgname.GO, 1, nchar(pkgname.GO) -3), "ENSEMBL", sep = ""))))) > length(unique(gwas$geneid)) / 2)
      geneid.type <- "ensembl"
  } else {
    if(sum(unique(gwas$genename) %in%  unlist(AnnotationDbi::as.list(get(paste(substr(pkgname.GO, 1, nchar(pkgname.GO) -3), "SYMBOL", sep = ""))))) > length(unique(gwas$geneid)) / 2)
      geneid.type <- "symbol"
    else if(sum(unique(gwas$genename) %in%  unlist(AnnotationDbi::as.list(get(paste(substr(pkgname.GO, 1, nchar(pkgname.GO) -3), "GENENAME", sep = ""))))) > length(unique(gwas$geneid)) / 2)
      geneid.type <- "genename"
  }
  if(is.null(geneid.type))
    stop("Cannot auto-detect gene identifier used (less than half of the GWAS gene identifies match the list of known gene identifiers for entrez gene IDs, ENSEMBL gene ids, symbols or genenames).")
  else if(geneid.type %in% c("entrez", "ensembl"))
    geneid.col <- "geneid"
  else if(geneid.type %in% c("symbol", "genename"))
    geneid.col <- "genename"
  message(geneid.type)
  
  # make gene identifiers unique, with lowest p being representative
  gwas <- gwas[order(gwas$gene.p), c(geneid.col, "gene.p")]
  gwas <- gwas[!duplicated(gwas[, geneid.col]), ]
  gwas <- na.omit(gwas)
  if(nrow(gwas) < 1) stop("Argument 'gwas' contains no data (or only NA rows)")
  
  allg <- gwas$gene.p
  names(allg) <- gwas[, geneid.col]
  
  message("Building topGO data object...")  
  capture.output(
      GO2genes <- annFUN.org(
          whichOnto = ontology, 
          mapping = pkgname.GO, 
          ID = geneid.type
      )
  )
  capture.output(
      tg.obj <- new(
          "topGOdata",
          ontology = ontology,
          allGenes = allg, 
          # in fact, geneSel is not important when using ks statistic (only used in the graph when showing gene counts per term, but we will have useInfo = 'np') 
          geneSel = function(foo) {foo < 0.01},
          nodeSize = pruneTermsBySize, 
          annot = annFUN.GO2genes,
          GO2genes = GO2genes
      )
  )
  print(tg.obj)
  
  message("Performing GO enrichment (delegated to the topGO package, using kolmogorov-smirnov test)...")
  capture.output(
      testdat <- runTest(
          tg.obj, 
          algorithm = topGOalgorithm, 
          statistic = "ks"
      )
  )
  print(testdat)
  
  if(!is.null(plotSigTermsToFile) && is.numeric(plotSigTermsToFile) && plotSigTermsToFile > 0) {
    message("\nPlotting graph (delegated to the topGO package)...\n")
    if(plotSigTermsToFile <= 0) {
      warning("Cannot plot enrichment graph - no significant terms available.")
    } else {
      if(geneData(testdat)["SigTerms"] <= plotSigTermsToFile) {
        plotSigTermsToFile <- geneData(testdat)["SigTerms"]
        warning(paste("Only", geneData(testdat)["SigTerms"], "terms available for plotting the enrichment graph."))
      }
      capture.output(suppressPackageStartupMessages(printGraph(tg.obj, testdat, firstSigNodes = plotSigTermsToFile, fn.prefix = "enrichment", pdfSW = TRUE, useInfo = c("np"))))
    }
  }
  
  # this is a numeric vector of p-values, named by GO term IDs
  res <- score(testdat)
  # annotate term names
  goids <- toTable(GOTERM)
  goids <- goids[!duplicated(goids$go_id), c("go_id", "Term")]
  res   <- merge(data.frame(go_id = names(res), P = res), goids, all.x = TRUE)
  res   <- vectorElements(res[order(res$P), ])
    
  return(res)
}
