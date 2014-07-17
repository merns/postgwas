# g is a graph
# value: list(graph = the colorized graph, legend.v = a data frame of vertices for the legend)
gwas2network.vertexcolor.overrep <- function(
  g, # igraph graph, vertex attribs 'name' and 'P'
  univ, # data frame of gene ids / names according to geneid.col and column P
  geneid.col, # either "geneid" or "genename" 
  pkgname.GO = "org.Hs.eg.db", 
  statistic = "fisher", # can also be "ks" for enrichment, see corresponding topGO argument 
  toFile = NULL
) {
  
  # set base color (has to be set because recolorization depends on the previous color)
  V(g)$color <- "black"
  signif.terms <- NULL
  legend.v <- NULL
  
  g.df <- data.frame(V(g)$name, V(g)$P)
  colnames(g.df) <- c(geneid.col, "P")
  
  tryCatch({
    signif.terms <- gwas2networkGOoverrep(
        univ = univ, 
        gwas = g.df, 
        geneid.col = geneid.col,
        pkgname.GO = pkgname.GO,
        # test for mult.testing many GO terms? Unsure because n is unkown especially for elim algorithm
        corr = FALSE, 
        statistic = statistic, 
        toFile = toFile
    )}, 
    error = function(e) {
      message(paste("ERROR: Could not set vertex color by GO overrepresentation analysis:  ", e))
    }
  )
  
  numterms <- if(is.null(signif.terms)) 0 else length(signif.terms$go_id)
  
  if(numterms <= 0) {
    message("No GO-terms with significant overrepresentation of GWAS genes found.")
  } else {
#        message("\nChoose terms to colorize? [Y/N]"); 
#        ans <- scan(n = 1, what = character(), quiet = TRUE); 
#        length(ans) > 0 && ans %in% c("Y", "y")
    
    # list: names = color, values = genes for that term
    clr.vnames <- signif.terms$genes
    names(clr.vnames) <- rainbow(numterms)

    if(statistic == "ks")
      signif.terms[["Plabel"]] <- paste("p (kolm. smirn.):", format.pval(signif.terms[["P"]]))
    else
      signif.terms[["Plabel"]] <- paste("p (fisher):", format.pval(signif.terms[["P"]]))
    
    # define legend vertices
    legend.v <- data.frame(
        name = paste("postgwasEnrichLegendVertex", c(paste(signif.terms[["go_id"]], signif.terms[["P"]], signif.terms[["Term"]], sep = "\n"), "brown: multiple overrepresented terms")),
        label = c(paste(signif.terms[["go_id"]], signif.terms[["Plabel"]], splitSentence(signif.terms[["Term"]]), sep = "\n"), "brown: multiple overrepresented terms"),
        # the weight is the overrepresentation p-value, multi-anotation is added with the minimum p-value of 1
        weight = c(as.numeric(signif.terms[["P"]]), 1), 
        P = NA, 
        SNP = NA, 
        degree = 0,
        marked = FALSE, 
        label.color = "black",
        color = c(names(clr.vnames), "saddlebrown")
    )
    
    for(clr in names(clr.vnames)) {
      V(g)$color[V(g)$color != "black" & V(g)$name %in% clr.vnames[[clr]]] <- "saddlebrown"
      V(g)$color[V(g)$color == "black" & V(g)$name %in% clr.vnames[[clr]]] <- clr
    }
      
  } # end signif terms found
  
  return(list(graph = g, legend.v = legend.v))
}




# value: new colorized graph (v$color set)
gwas2network.vertexcolor.regex <- function(
  g, # igraph graph with vertex attribute name
  regex, 
  geneid.col, # name of the column that is used in the graph for names 
  bm.config = NULL, 
  ds = bm.init.genes(bm.config),
  use.buffer = FALSE
  ) {
  
  if(!all(names(regex) %in% colors()))
    stop("Vertex color regex list contains invalid color names. See colors() for all valid color names.\n")
  
  message("Retrieving gene ontology annotation for vertex color...")
  go <- bm.genes.goterms(
      config = bm.config, 
      ds = ds, 
      filter.name = if(geneid.col == "genename") bm.config$gene$filter$name else bm.config$gene$filter$id, 
      filter.val = V(g)$name,
      use.buffer = use.buffer
  )
  
  # names : color, value: a list of genenames (vertices) to colorize
  clr.vnames <- lapply(
    regex, 
    function(expr) unique(go[grep(expr, as.vector(go$"goterm.name"), ignore.case = T), geneid.col])
    )
  
  V(g)$color <- "black"
  
  for(clr in names(clr.vnames)) {
    V(g)$color[V(g)$color != "black" & V(g)$name %in% clr.vnames[[clr]]] <- "saddlebrown"
    V(g)$color[V(g)$color == "black" & V(g)$name %in% clr.vnames[[clr]]] <- clr
  }
  
  # define legend vertices
  legend.v <- data.frame(
      name = paste("postgwasEnrichLegendVertex", c(paste("GO term name matches: \"", regex, "\"", sep = ""), "brown: multiple\nmatching terms")),
      label = c(splitSentence(paste("GO term name matches: \"", regex, "\"", sep = "")), "brown: multiple\nmatching terms"),
      weight = 1, 
      P = NA, 
      SNP = NA, 
      degree = 0,
      marked = FALSE, 
      label.color = "black",
      color = c(names(regex), "saddlebrown")
  )
  
  return(list(graph = g, legend.v = legend.v))
}
