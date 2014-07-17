# returns NULL when overrep analysis does not yield (valid) results
gwas2networkGOoverrep <- function(
    univ, # data frame of geneids or symbols, according to geneid.col, and a P column 
    gwas, # data frame with columns: P and value of geneid.col 
    geneid.col = "geneid",
    pkgname.GO = "org.Hs.eg.db",
    corr = TRUE,     # test for mult.testing many GO terms? Unsure because n is unkown especially for elim algorithm
    statistic = "fisher", # can also be "ks" for enrichment instead of overrepresentation, see corresponding topGO argument 
    toFile = "gwas2network"
) {
  
  # prepare gene sets:
  # make gene identifiers unique, with lowest p being representative
  gwas <- gwas[order(gwas$P), c(geneid.col, "P")]
  gwas <- gwas[!duplicated(gwas[, geneid.col]), ]
  # remove results where either identifier or P is NA
  gwas <- na.omit(gwas)
  if(nrow(gwas) < 1)
    return(NULL)
  # need only gene identifier for group
  grp <- gwas[, geneid.col]
  
  # similar procedure for gene universe
  univ <- na.omit(univ)
  # make univ genes unique, with lowest p being representative
  univ <- univ[order(univ$P), ]
  univ <- univ[!duplicated(univ[, geneid.col]), ]
  # finally, create named vector of p-values
  allg <- univ$P
  names(allg) <- univ[, geneid.col]
  
  # build topGO data object
  message("Building topGO data object...")
  capture.output(
      GO2genes <- annFUN.org(
          whichOnto = "BP", 
          mapping = pkgname.GO, 
          ID = if(geneid.col == "geneid") "entrez" else "symbol"
      )
  )
  if(sum(names(allg) %in% unlist(GO2genes)) < length(allg) / 2)
    stop("More than half of genes can not be annotated to GO terms. Probably gene identifier / names from the arguments are not compatible with the GO annotation package used. Normally, using entrez IDs (in the column 'geneid') or gene symbols (column 'genename') should work. ")
  capture.output(
      tg.obj <- new(
          "topGOdata",
          ontology = "BP",
          allGenes = allg, 
          geneSel = function(uselessScoreVec) {names(allg) %in% grp},
          annot = annFUN.GO2genes,
          GO2genes = GO2genes
      )
  )
  
  if(statistic == "ks")
    message("Performing GO enrichment (topGO)...")
  else
    message("Performing GO overrepresentation analysis (topGO)...")
  capture.output(
      testdat <- runTest(
          tg.obj, 
          algorithm = "weight01", 
          statistic = statistic
      )
  )
  
  if(!is.null(toFile) && is.character(toFile)) {
    capture.output(printGraph(tg.obj, testdat, firstSigNodes = 3, fn.prefix = toFile, pdfSW = TRUE, useInfo = "all"))
  }
  
  if(geneData(testdat)["Significant"] > 0 && geneData(testdat)["SigTerms"] > 0) {
    res <- na.omit(sort(score(testdat), decreasing = FALSE)[1:3])
    # get mapping of GO IDs to name etc as data frame from the org... package
    goids <- toTable(GOTERM)
    goids <- goids[!duplicated(goids$go_id), c(2, which(names(goids) == "Term"))]
    # annotate term names to results
    res <- merge(goids, data.frame(go_id = names(res), P = res), all.y = TRUE)
    res <- vectorElements(res[order(res$P), ])
    
    return(list(
            go_id = res$go_id, 
            # if correction for multiple tests is requested, length(usedGO(tg.obj)) is the number of GO terms used / tested
            # this does not account for the elim algorithm which tests probably only a subset of these (=> p is conservative / less power)
            P = if(corr) p.adjust(res$P, n = length(usedGO(tg.obj))) else res$P,
            Term = res$Term,
            # all terms from res are in GO2genes, because re has been computed using GO2genes
            genes = GO2genes[res$go_id]
        ))
  } else {
    return(NULL)
  }
}







# test for an accumulation of high-evidence association scores in functionally connected modules 
gwas2networkModuleEnrich <- function(
    wholegraph, # igraph graph object with vertex attributes name and P
    module # igraph graph object with vertex attributes name and P
) {
  
  whole <- V(wholegraph)$P
  names(whole) <- V(wholegraph)$name
  
  mod <- V(module)$P
  names(mod) <- V(module)$name
  
  # all genes in the module have to be present in the universe network
  if(!all(names(mod) %in% names(whole))) {
    warning("gwas2networkModuleEnrich; Genes from the module are not contained in the global network. Enrichment score is not being calculated.")
    p <- NA
  } else {
    dat <- data.frame(p = V(wholegraph)$P, inmod = V(wholegraph)$name %in% V(module)$name)
    p <- wilcox.test(p ~ inmod, data = dat, alternative = "greater")$p.value
    # greater/less relates to the first argument, we are interested whether mod p-values are smaller (less) 
#    p <- wilcox.test(mod, whole, alternative = "less")$p.value
  }
  
  legend.v <- data.frame(
      name = "postgwasEnrichLegendVertexModule", 
      label = paste("Module score:", signif(p)),
      # has lowest weight and will be listed first in the legend
      weight = 0, 
      P = NA,
      SNP = NA, 
      degree = 0,
      marked = FALSE, 
      label.color = "black",
      color = "black"
  ) 
  
  return(list(p = p, legend.v = legend.v)) 
}
