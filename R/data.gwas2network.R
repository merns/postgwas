utils::globalVariables(c("gwas2network.genename.x", "gwas2network.genename.y", "weight", ".SD"))

# assert that at execution time, the 'network' varaible has been set in the postgwasBuffer environment with current network data.
# prepares the vertex and edge data and returns it
# defines vertex columns marked, gwas, degree, weight, label, label.color and category columns
# defines edge columns weight, color
# marked column indicate dual annotations (single association p-value given to several genes)
data.gwas2network <- function(
  gwas.mapped.genes, 
  geneid.col, 
  prune, 
  remove.superhubs,
  p.default, 
  edge.weight.fun,
  edge.weight.collapse.fun,
  toFile = TRUE
) {

  gwas.mapped.genes[, geneid.col] <- as.vector(toupper(gwas.mapped.genes[, geneid.col]))
  genes.gwas <- unique(as.vector(gwas.mapped.genes[, geneid.col]))
  
  # get a local copy of the (readily processed) network data
  e <- get("network", envir = postgwasBuffer)
  
  if(!is.null(prune)) {
    message(paste("Reducing network to", prune, "interactors"))
    if(prune == "gwasonly") {
      e <- e[e[, 1] %in% genes.gwas & e[, 2] %in% genes.gwas, ]
    }
    if(prune %in% c("connected", "shared")) {
      e <- e[e[, 1] %in% genes.gwas | e[, 2] %in% genes.gwas, ]
    }
    # for shared pruning, edges have to be unique - we do that later
    if(nrow(e) < 1)
      stop("Network has no edges after pruning!\n")
  }
  
  # make edges unique
  message("Searching for duplicate edges (if existing) and collapsing these to single edges: ")
  colnames(e)[1:2] <- c("gwas2network.genename.x", "gwas2network.genename.y")
  e$weight <- as.numeric(as.vector(e$weight))
  e <- data.table(e)
  setkey(e, gwas2network.genename.x, gwas2network.genename.y)
  
  message("Collapsing weights (mean by default) ... ")
  e <- e[TRUE, weight := edge.weight.collapse.fun(.SD$weight, .SD$label), by = c("gwas2network.genename.x", "gwas2network.genename.y")]

  message("Collapsing labels ... ")
  if(!is.null(e$label)) { 
    e <- e[
        TRUE, 
        label := {lbls <- unique(.SD$label); if(length(lbls) > 3) "> 3 labels" else paste(lbls, collapse = " / ")}, 
        by = c("gwas2network.genename.x", "gwas2network.genename.y")
    ]
  }
  
  # keep only first of all duplicate edges
  e <- as.data.frame(e[!duplicated(e), ])

  # remove loops
  e <- unique(e[e[, 1] != e[, 2], ])
  # also remove empty string identifier - causes problems
  e <- e[e[, 1] != "" & e[, 1] != "", ]
  
  if(!is.null(prune) && prune == "shared") {
    # remove all edges with genes that have only one interactor and are not gwas genes
    ecounts <- table(as.vector(unlist(e[, 1:2])))
    remove <- names(ecounts[ecounts <= 1 & !(names(ecounts) %in% genes.gwas)])
    e <- e[!(e[, 1] %in% remove | e[, 2] %in% remove), ]
  }
  
  # normally, it is a good idea to remove super-hubs (v > mean(degree)^2)
  # we do this after pruning to allow valid auto-detection with network size (default argument remove.superhubs)
  if(remove.superhubs && prune != "gwasonly") {
    v <- data.frame(name = unique(as.vector(unlist(e[, 1:2]))), stringsAsFactors = FALSE)
    testg <- graph.data.frame(e, directed = FALSE, v)
    superhubs <- V(testg)$name[igraph::degree(testg) > mean(igraph::degree(testg))^2]
    superhubs <- superhubs[!(superhubs %in% genes.gwas)]
    message(paste("Removing superhubs:", paste(superhubs, collapse = " ; ")))
    e <- e[!(e[, 1] %in% superhubs), ]
    e <- e[!(e[, 2] %in% superhubs), ]
  }
  
  if(nrow(e) < 1)
    stop("Network has no edges after pruning and processing!\n")
  else
    message("Network has ", nrow(e), " edges")
  

  ######### vertex weight, label and edge attributes #########
  
  message("Setting vertex weights...")
  # order genes by p-value and remove duplicated genes (will keep the smallest SNP p-value for each gene)
  # vertices are all in network and all in gwas
  v <- vectorElements(merge(data.frame(name = unique(c(genes.gwas, as.vector(unlist(e[, 1:2]))))), gwas.mapped.genes, by.x = "name", by.y = geneid.col, all.x = T))
  v$P <- as.numeric(v$P)
  if(is.null(v$pheno))
    v$pheno <- "defaultPhenotype"
  v$pheno[is.na(v$pheno)] <- "defaultPhenotype"
  
  # distinguish by pheno - for each pheno the smallest P per gene
  # make genes for each pheno unique with lowest P
  v <- list2df(by(
    v, 
    v$pheno,
    function(v.phen) {
      v.phen <- v.phen[order(v.phen$P), ]
      v.phen <- v.phen[!duplicated(v.phen$name), ]
      # check for dual annot within this pheno (vertices that got their weight from the same SNP)
      v.phen$marked <- v.phen$SNP %in% unique(as.vector(v.phen[duplicated(v.phen$SNP) & !is.na(v.phen$SNP), "SNP"]))
      # write dual annot file
      marked <- gwas.mapped.genes[gwas.mapped.genes$SNP %in% v.phen[v.phen$marked, "SNP"], ]
      if(toFile && nrow(marked) > 0) {
        marked <- by(marked[, geneid.col], factor(marked$SNP), function(x) paste(unique(x), collapse = ";"))
        writeLines(paste(names(marked), marked), paste("gwas2networkMultiAnnotations", unique(v.phen$pheno), ".txt", sep = ""))
      }
      return(v.phen)
    },
    simplify = FALSE
  ))
  
  v[is.na(v$P), "P"] <- p.default
  v$weight <- -log10(v$P)

  # make genes globally unique by collapsing genes with multiple phenos in additional columns
  # the lowest P will be the representative P (in the P column), further phenotypes are shifted to the additional columns
  # we rename phenotypes with predefined plot symbols
  phenos <- names(table(v$pheno)[order(table(v$pheno), decreasing = T)]) # unique phenos, ordered by abundance
  if(length(phenos) > 3)
    stop("Too many phenotypes - cannot plot that (note that vertices in the network that are not in gwas.mapped.genes form a dummy phenotype).\n")
  symbols <- c("circle", "square", "diamond")
  for(i in 1:length(phenos))
    v$pheno[v$pheno == phenos[i]] <- symbols[i]
  # set additional columns with defaults
  v[, c(unique(v$pheno), paste(unique(v$pheno), "marked", sep = "."), paste(unique(v$pheno), "weight", sep = "."))] <- FALSE

  v <- list2df(by(
    v, 
    factor(v$name),
    function(v.name) {
      v.name <- v.name[order(v.name$P), ]
      single <- v.name[1, ] # lowest P is default pheno
      if(nrow(v.name) > 1) {
        # when one pheno has dual annot, set true
        # this is mandatory because we later use the set of dual-annot genes
        single$marked <- sum(v.name$marked) > 0
        # now set collapsing columns for additional rows
        for(i in 2:nrow(v.name)) {
          curr.phen <- v.name$pheno[i]
          single[1, curr.phen] <- TRUE
          single[1, paste(curr.phen, "marked", sep = ".")] <- v.name$marked[i]
          single[1, paste(curr.phen, "weight", sep = ".")] <- v.name$weight[i]
        }
      }
      return(single)
    }, 
    simplify = FALSE
  ))
  
  # phenotypes -> categories
  colnames(v)[colnames(v) == "pheno"] <- "category.main"
  
  v$degree <- igraph::degree(graph.data.frame(e, directed = FALSE, v))
  v$gwas   <- v$name %in% genes.gwas
  
  v$label  <- v$name
  
  
  ######### edge weights #########
  
  message("Setting edge weights...")
  colnames(e)[1:2] <- c("gwas2network.genename.x", "gwas2network.genename.y")
  colnames(e)[colnames(e) == "weight"] <- "net.weight.temp"
  # assign vertex p-values to each edge
  ew <- merge(e, v, by.x = "gwas2network.genename.x", by.y = "name", all.x = T)
  ew <- merge(ew, v, by.x = "gwas2network.genename.y", by.y = "name", all.x = T)
  # now we have in ew the columns degree.x, degree.y and P.x and P.y
  e <- e[order(e$gwas2network.genename.x, e$gwas2network.genename.y), ]
  ew <- ew[order(ew$gwas2network.genename.x, ew$gwas2network.genename.y), ]
  e$weight <- edge.weight.fun(p1 = ew$P.x, p2 = ew$P.y, degree1 = ew$degree.x, degree2 = ew$degree.y, weight = as.numeric(ew$net.weight.temp))
  
  
  ######### write parameters #########
  
  if(toFile) {
    writeLines(
      c(
        paste("Default P =", p.default), 
        paste("gwas.mapped.genes had", nrow(gwas.mapped.genes), "rows"),
        paste("edge weight function for community detection: ", paste(deparse(body(edge.weight.fun)), sep = "\n")),
        if(remove.superhubs && prune != "gwasonly") paste("Superhubs removed:", paste(superhubs, collapse = " ; "))
        ),
      "gwas2networkParams.txt"
      )
  }
  
  return(list(e = e, v = v))
}
