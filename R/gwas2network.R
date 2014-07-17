gwas2network <- function(
    gwas.mapped.genes, 
    network, 
    prune = "gwasonly",
    max.communities = 500, 
    vertexcolor.GO.overrep = "org.Hs.eg.db",
    vertexcolor.GO.regex = NULL,
    min.transparency = 0.25, 
    max.transparency = 0.75, 
    custom.layout = FALSE,
    edge.weight.fun = function(p1, p2, degree1, degree2, weight) {
      return(((-log10(p1) +1) * (-log10(p2) +1))^2 * weight)
    },
    edge.weight.collapse.fun = function(weights, labels) {
      return(mean(weights))
    },
    remove.superhubs = nrow(network) > 100,
    p.default = max(gwas.mapped.genes$P, na.rm = TRUE),
    png.maxedges = 2500, 
    file.verbosity = 2, 
    biomart.config = biomartConfigs$hsapiens, 
    use.buffer = FALSE, 
    cores = 1
){
  
  if(is.null(biomart.config) || !is.list(biomart.config))
    stop("Argument 'biomart.config' has to be a list")
  
  if(!is.function(edge.weight.fun) || !is.function(edge.weight.collapse.fun))
    stop("Arguments 'edge.weight.fun' and 'edge.weight.collapse.fun' have to be functions.\n")
  
  if(!all(c("p1", "p2", "degree1", "degree2", "weight") %in% formalArgs(edge.weight.fun)))
    stop("Argument 'edge.weight.fun' has an incomplete argument list, must contain 'p1', 'p2', 'degree1', 'degree2', 'weight'.\n")
  
  if(!all(c("weights", "labels") %in% formalArgs(edge.weight.collapse.fun)))
    stop("Argument 'edge.weight.collapse.fun' has an incomplete argument list, must contain 'weights', 'labels'.\n")
  
  if(!is.null(cores) && cores > 1) {
    if(!'parallel' %in% installed.packages()[, 'Package'])
      stop("Function parameter cores > 1 requires the package 'parallel' to work correctly.\n")
    suppressPackageStartupMessages(stopifnot(require(parallel, quietly = TRUE)))
    myapply <- mclapply
  } else {
    myapply <- lapply
  }   
  
  if(!is.null(vertexcolor.GO.regex)) {
    if(!is.list(vertexcolor.GO.regex))
      stop("Argument 'vertexcolor.GO.regex' is not a list.\n")
    if(!all(names(vertexcolor.GO.regex) %in% colors()))
      stop("Argument 'vertexcolor.GO.regex' contains invalid color names. See colors() for a list of valid color names.\n")
  }
  
  if(is.null(vertexcolor.GO.regex) && !is.null(vertexcolor.GO.overrep)) {
    if(!is.vector(vertexcolor.GO.overrep) || length(vertexcolor.GO.overrep) > 1 || typeof(vertexcolor.GO.overrep) != "character")
      stop("Argument 'vertexcolor.GO.overrep' has to be a character string defining an Annotation package name (see help pages).\n")
    if(!("topGO" %in% installed.packages()[, "Package"])) {
      if(interactive() && readline("Package 'topGO' is not installed. Try to install? [Y/N] ") %in% c("Y", "y")) {
        biocLite <- NULL
        source("http://bioconductor.org/biocLite.R")
        biocLite("topGO")
      } else {
        stop("Vertex colorization by GO overrepresentation analsis requires the topGO package to be installed.\n")		
      }
    }
    if(biomart.config$gene$attr$id != "entrezgene") 
      warning("Vertex colorization by significance normally requires entrez gene ids (specify properly in the biomart configuration)\n")
    require(topGO)
    # load proper GO organism data
    success <- do.call(library, list(package = vertexcolor.GO.overrep, logical.return = TRUE))
    if(!success)
      stop("The GeneOnotology annotation package supplied in argument 'vertexcolor.GO.overrep' is not installed or cannot be loaded.\n")
  }
  
  # there is a legacy package 'igraph0' that masks igraph and leads to errors
  # ensure that our igraph is first on search path 
  # and remember its current position for restore at function exit
#  igraph.pos.orig <- which(search() == "package:igraph")
#  suppressWarnings(detach("package:igraph", force = TRUE))
#  suppressWarnings(suppressPackageStartupMessages(library(igraph, pos = 2, quietly = TRUE)))
  
  
  if(is.null(network) || nrow(network) == 0)
    stop("Argument 'network' does not contain data.\n")
  
  if(ncol(network) < 2)
    stop("Argument 'network' has to contain at least two columns with vertex IDs.\n")
  if(any(colnames(network)[1:2] %in% c("weight", "label")))
    stop("The first or second column of the 'network' argument is named 'weight' or 'label'.\n")
  
  if(any(c("net.weight.temp", "gwas2network.genename.x", "gwas2network.genename.y") %in% colnames(network))) 
    stop("Parameter 'network' may not contain the columns 'net.weight.temp', 'gwas2network.genename.x'', 'gwas2network.genename.y'.\n")
  
  network.allgenes <- unique(c(as.vector(network[, 1]), as.vector(network[, 2])))
  
  # because the 'network' argument is potentially large, we put it in the postgwasBuffer enviroment
  # so that is does not have to be passed as function argument
  assign("network", network, envir = postgwasBuffer)
  remove.superhubs # we have to evaluate the remove.superhubs argument (that depends on 'network') before the 'network' variable vanishes
  rm(network)
  
  if(is.null(gwas.mapped.genes) || !is.data.frame(gwas.mapped.genes) || nrow(gwas.mapped.genes) == 0) {
    stop("Argument 'gwas.mapped.genes' has to be a nonempty data frame.\n")
  } else {
    gwas.mapped.genes <- vectorElements(gwas.mapped.genes)
    if(any(c("net.weight.temp", "gwas2network.genename.x", "gwas2network.genename.y") %in% colnames(gwas.mapped.genes))) 
      stop("Argument 'gwas.mapped.genes' may not contain the columns 'net.weight.temp', 'gwas2network.genename.x'', 'gwas2network.genename.y'.\n")
    if(!is.null(gwas.mapped.genes$gene.p)) {
      message("Using gene-based p-values (column gene.p)")
      gwas.mapped.genes$P <- gwas.mapped.genes$gene.p
    } else {
      if(is.null(gwas.mapped.genes$P)) 
        stop("Argument 'gwas.mapped.genes' has to contain a column 'P' or 'gene.p'.\n")
    }
    if(!any(c("geneid", "genename") %in% colnames(gwas.mapped.genes))) 
      stop("Argument 'gwas.mapped.genes' has to contain either a column 'geneid' or 'genename'.\n")
    if(all(c("geneid", "genename") %in% colnames(gwas.mapped.genes))) {
      # determine the kind of identifier used in the network data
      if(sum(as.vector(na.omit(gwas.mapped.genes[, "geneid"])) %in% network.allgenes)
          > 
         sum(as.vector(na.omit(gwas.mapped.genes[, "genename"])) %in% network.allgenes)
         )
        geneid.col <- "geneid"
      else
        geneid.col <- "genename"
    } else {
      geneid.col <- if("geneid" %in% colnames(gwas.mapped.genes)) "geneid" else "genename"
    }
    message(paste(
      "Using", 
      geneid.col, 
      "as primary identifier.", 
      round(sum(unique(gwas.mapped.genes[, geneid.col]) %in% network.allgenes) / length(unique(gwas.mapped.genes[, geneid.col])) * 100, 2),
      "% of query genes are contained in the network before processing."
    ))
  }
  
  # format the network data 
  data.gwas2network.preprocessNet()
  
  g.dat <- data.gwas2network(
              gwas.mapped.genes = gwas.mapped.genes, 
              geneid.col = geneid.col,
              prune = prune, 
              remove.superhubs = remove.superhubs,
              p.default = p.default, 
              edge.weight.fun = edge.weight.fun,
              edge.weight.collapse.fun = edge.weight.collapse.fun,
              toFile = file.verbosity > 1
            )
  e <- g.dat$e
  v <- g.dat$v
  rm(g.dat)
  
  # we need an idmap for two purposes:
  # when ids are used in the graph, we want to display names in the plot
  # we need a universe of all genes for the organism, which is represented by the (complete) idmap
  message("Retrieving all gene IDs and names from biomart for the organism...")
  delayedAssign("ds", bm.init.genes(biomart.config))
  idmap <- bm.genes(
      config = biomart.config, 
      ds = ds, 
      use.buffer = use.buffer
  )[, 1:2]
  idmap.nona <- na.omit(idmap)
  # add P when known (discards genes that are not in univ)
  idmap.p <- merge(idmap.nona, data.frame(name = v$name, P = v$P), all.x = TRUE, by.x = geneid.col, by.y = "name")
  # unknown p's to default
  idmap.p[is.na(idmap.p$P), "P"] <- p.default

  
  g <- graph.data.frame(e, directed = FALSE, v)
  

  ######### vertex color and graph plot #########
  
  if(!is.null(vertexcolor.GO.overrep) || !is.null(vertexcolor.GO.regex)) {
    message("Setting vertex color...")
    if(!is.null(vertexcolor.GO.regex)) {
      clr.res <- gwas2network.vertexcolor.regex(
        g = g, 
        regex = vertexcolor.GO.regex, 
        geneid.col = geneid.col,
        bm.config = biomart.config,
        ds = ds, 
        use.buffer = use.buffer
      )
    } else {
      clr.res <- gwas2network.vertexcolor.overrep(
        g = g, 
        univ = idmap.p,
        geneid.col = geneid.col, 
        pkgname.GO = vertexcolor.GO.overrep,
        toFile = if(file.verbosity > 1) "gwas2networkGOoverrep" else NULL
      )
    }
    g <- clr.res$graph
    legend.v <- clr.res$legend.v
  } else {
    # no GO colors
    V(g)$color <- "black"
    legend.v <- NULL
  }
  
  
  ######### set transparency #########

  vw.max <- max(V(g)$weight)
  vw.min <- min(V(g)$weight)
  ew.max <- max(E(g)$weight)
  ew.min <- min(E(g)$weight)

  if(vw.min < 0)
    stop("Vertex weights < 0 found. Aborting.\n")

  # determine alpha (normalize weight to [0,1] and then expand to [min.transparence,max.transparency])
  if(vw.max == vw.min) {
    # avoid division by NULL
    V(g)$alpha <- max.transparency
  } else {
    V(g)$alpha <- ((V(g)$weight - vw.min) / (vw.max - vw.min)) * (max.transparency - min.transparency) + min.transparency
  }
  
  if(ew.max == ew.min) {
    E(g)$alpha <- max.transparency
  } else {
    E(g)$alpha <- ((E(g)$weight - ew.min) / (ew.max - ew.min)) * (max.transparency - min.transparency) + min.transparency
  }
  
  V(g)$label.color <- rgb(red = 0, green = 0, blue = 0, alpha = V(g)$alpha)
  E(g)$color       <- rgb(red = 0, green = 0, blue = 0, alpha = E(g)$alpha^1.7)
  # the vertex transparency has to be adjusted for each graph seperately, because the base vertex color can change after GO overrep analysis
  
  
  ######### define scale vertices #########
  # scale vertices do not have valid gene names and may only be added after GO overrep
  scale.vertices <- list(nv = 0)
  if(vw.max != vw.min) {
    bins <-scaleBins(min(vw.min, 2), max(2, ceiling(vw.max)), 5, 5)
    alpha <- ((bins -vw.min) / (vw.max - vw.min)) * (max.transparency - min.transparency) + min.transparency
    alpha[alpha > 1] <- 1
    alpha[alpha < 0] <- 0
    color <- rgb(0, 0, 0, alpha = alpha)
    scale.vertices <- list(
      nv = length(bins),
      name = paste("postgwasScaleVertex", bins, sep = ""), 
      label = format.pval(10^(-(bins))),
      weight = bins, 
      P = 10^(-(bins)), 
      SNP = format.pval(10^(-(bins))),
      marked = FALSE,
      label.color = color,
      color = color
    )
  }
  
  
  message("Plotting the entire graph...")
  g.whole <- g
  v.color          <- t(col2rgb(V(g.whole)$color)) / 255
  V(g.whole)$color <- rgb(red = v.color[, "red"], green = v.color[, "green"], blue = v.color[, "blue"], alpha = V(g.whole)$alpha)
  # translate labels to names when IDs were used
  if(geneid.col == "geneid")
    g.whole <- graph.id2name(g.whole, idmap.nona)
  g.whole <- addLegendVertices(g.whole, legend.v)
  # add scale vertices
  scale.vertices$graph <- g.whole
  g.whole <- do.call(add.vertices, scale.vertices)
  if(file.verbosity > 0) {
    gwas2network.plot(g.whole, file = "gwas2network.pdf", title = "gwas2network plot", v.weight.max = vw.max)
    if(ecount(g.whole) <= png.maxedges)
      gwas2network.plot(g.whole, device = png, file = "gwas2network.png", res = 200, v.weight.max = vw.max)
  }

  
  ######### community detection and enrichment (GO/p-values) #########
  
  g.comps <- decompose.graph(g, mode = "strong", min.vertices = 2)  
  g.comps <- g.comps[order(sapply(g.comps, vcount), decreasing = TRUE)]
  
  if(max.communities > 1) {
    message("Detecting communities and calculating community scores ...")
    subgs <- list()
    # detect communities in each connected component seperately and collect them in subgs
    for(g.comp in g.comps) {
      membership <- spinglass.community(
          graph = g.comp, 
          weights=E(g.comp)$weight, 
          spins = max.communities,
          implementation = ifelse(ew.min < 0, "neg", "orig")
      )$membership
      subgs.comp <- myapply(
          names(table(factor(membership))),
          function(mod.id, mc.cores) { # mc.cores is a dummy argument to catch the argument when normal apply is used
            subg <- induced.subgraph(g.comp, V(g.comp)[membership == mod.id])
            # discard empty and single-vertex communities which spinglass can produce both
            if(vcount(subg) <= 1)
              return(NA)
            # calculate p-value enrichment of community
            mod.enr <- gwas2networkModuleEnrich(g, subg)
            return(list(subg = subg, legend.v = mod.enr$legend.v, p = mod.enr$p))
          },
          mc.cores =cores
      )
      subgs.comp <- subgs.comp[!is.na(subgs.comp)]
      # collect for all connected components
      subgs <- c(subgs, subgs.comp)
    }
    
    # sort communities by p-value (for plotting and file names in correct order)
    subgs <- subgs[order(as.numeric(list2df(subgs)[, "p"]))]
    
    # vertexcolor by community specific GO term overrepresentation
    if(is.null(vertexcolor.GO.regex) && !is.null(vertexcolor.GO.overrep)) {
      subgs <- myapply(
          1:length(subgs), 
          function(idx, mc.cores) {
            subg <- subgs[[idx]]$subg
            legend.v <- subgs[[idx]]$legend.v
            V(subg)$color <- "black"
            message(paste("Setting vertex color of community", idx, "..."))
            # the universe with regard to the extracted modules is the entire graph (not all genes in the organism, because we sample, i.e. build communities from the graph)
#            univ <- data.frame(V(g)$name, V(g)$P)
#            colnames(univ) <- c(geneid.col, "P")
            clr.res <- gwas2network.vertexcolor.overrep(
                g = subg,
                univ = idmap.p, 
                geneid.col = geneid.col, 
                pkgname.GO = vertexcolor.GO.overrep, 
                toFile = if(file.verbosity > 2) paste("gwas2networkGOoverrepCommunity", idx, sep = "") else NULL
            )
            legend.v <- rbind(legend.v, clr.res$legend.v)
            subg <- clr.res$graph
            # set transparency 
            v.color       <- t(col2rgb(V(subg)$color)) / 255
            V(subg)$color <- rgb(red = v.color[, "red"], green = v.color[, "green"], blue = v.color[, "blue"], alpha = V(subg)$alpha)
            return(list(subg = subg, legend.v = legend.v))
          },
          mc.cores =cores
      )
    }
    
    # add legend and scale to communities
    # also translate vertex labels to names (when IDs were used)
    subgs <- myapply(
        1:length(subgs), 
        function(idx, mc.cores) {
          subg <- subgs[[idx]]$subg
          if(geneid.col == "geneid")
            subg <- graph.id2name(subg, idmap.nona)
          subg <- addLegendVertices(subg, subgs[[idx]]$legend.v)
          scale.vertices$graph <- subg
          subg <- do.call(add.vertices, scale.vertices)
          return(subg)
        },
        mc.cores =cores
    )
    
    # finally, print and plot communities
    message("Plotting communities...")
    if(custom.layout && interactive() && cores > 1)
      mylayoutapply <- lapply # custom layout by drag and drop cannot be done parallel
    else
      mylayoutapply <- myapply 
    
    mylayoutapply(
        1:length(subgs), 
        function(idx, mc.cores) {
          subg <- subgs[[idx]]
          # dump subgraph data tables
          if(file.verbosity > 2) {
            v.dump <- merge(v, data.frame(name = V(subg)$name, color = V(subg)$color))
            write.table(v.dump, paste("gwas2networkGraphVerticesCommunity", idx, ".csv", sep = ""), row.names = FALSE, sep = "\t")
            e.dump <- e[e$gwas2network.genename.x %in% as.vector(get.edgelist(subg)) & e$gwas2network.genename.y %in% as.vector(get.edgelist(subg)), ]
            write.table(e.dump, paste("gwas2networkGraphEdgesCommunity", idx, ".csv", sep = ""), row.names = FALSE, sep = "\t")
          }
          # plot communities
          if(custom.layout && interactive() && ecount(subg) > 0) {
            message("Drag vertices to the custom layout. When done, press Enter (before closing the window).")
            subg.tkp <- subg
            V(subg.tkp)$frame.color <- V(subg.tkp)$color
            V(subg.tkp)$label.color <- V(subg.tkp)$color
            V(subg.tkp)$color <- "white"
            tkp.id <- suppressWarnings(tkplot(subg.tkp))
            scan(n = 1, quiet = T)
            layout <- tkplot.getcoords(tkp.id)
            tkplot.close(tkp.id)
            layout[, 2] <- -layout[, 2]
            layout <- layout / 2.5
          } else {
            layout <- layout.fruchterman.reingold(subg, area = vcount(subg)^3)
          }
          if(file.verbosity > 0) {
            gwas2network.plot(subg, layout = layout, title = "gwas2network plot (communities)", file = paste("gwas2networkCommunity", idx, ".pdf", sep = ""), v.weight.max = vw.max)
            if(ecount(subg) <= png.maxedges)
              gwas2network.plot(subg, layout = layout, device = png, file = paste("gwas2networkCommunity", idx, ".png", sep = ""), res = 200, v.weight.max = vw.max)
          }
        },
        mc.cores =cores
    )
    
  } else {
    # no community detection: each connected component is treated as community, without overrep
    message("Plotting connected components...")
    for(idx in 1:length(g.comps)) {
      g.comp <- g.comps[[idx]]
      # translate labels to names (when IDs were used)
      if(geneid.col == "geneid")
        g.comp <- graph.id2name(g.comp, idmap.nona)
      # transparency
      v.color         <- t(col2rgb(V(g.comp)$color)) / 255
      V(g.comp)$color <- rgb(red = v.color[, "red"], green = v.color[, "green"], blue = v.color[, "blue"], alpha = V(g.comp)$alpha)
      g.comp <- addLegendVertices(g.comp, legend.v)
      # add scale vertices
      scale.vertices$graph <- g.comp
      g.comp <- do.call(add.vertices, scale.vertices)
      if(file.verbosity > 0) {
        gwas2network.plot(g = g.comp, file = paste("gwas2networkConnectedComponent", idx, ".pdf", sep = ""), title = "gwas2network plot (connected components)", v.weight.max = vw.max)
        if(ecount(g.comp) <= png.maxedges) 
          gwas2network.plot(g = g.comp, file = paste("gwas2networkConnectedComponent", idx, ".png", sep = ""), device = png, res = 200, v.weight.max = vw.max)
      }
    }
  }
  
  
  ######### dump (whole) graph data, return #########

  v.dump <- merge(v, data.frame(name = V(g)$name, color = V(g)$color, size = v$weight * 2))
  if(file.verbosity > 1) {
    write.table(v.dump, "gwas2networkGraphVertices.csv", row.names = FALSE, sep = "\t")
    write.table(e, "gwas2networkGraphEdges.csv", row.names = FALSE, sep = "\t")
  }
  
  # restore search path
#  suppressWarnings(detach("package:igraph", force = TRUE))
#  suppressWarnings(suppressPackageStartupMessages(library(igraph, pos = igraph.pos.orig, quietly = TRUE)))
  
  return(g.whole)
  
}




addLegendVertices <- function(
                       graph, 
                       v.df, 
                       weight = mean(V(graph)$weight)
                     ) {
  if(is.null(v.df))
    return(graph)
  else
    return(add.vertices(
      graph, 
      nv = nrow(v.df),
      name = as.vector(v.df$name), 
      label = as.vector(v.df$label),
      weight = weight, 
      P = 10^(-weight), 
      SNP = v.df$label,
      marked = FALSE,
      label.color = as.vector(v.df$color),
      color = rgb(red = 0, green = 0, blue = 0, alpha = 0)
    ))
}

graph.id2name <- function(graph, idmap.nona) {
    V(graph)$label <- sapply(
                       V(graph)$label, 
                       function(x) {
                         mapped <- idmap.nona[idmap.nona$geneid == x, "genename"]
                         if(length(mapped) == 0 || mapped == "") x else mapped[1]
                       }
                     )
  return(graph)
}


# returns a reasonable number of integer bins (tick marks, ...) for an integer interval
# target bincounts can be equal, target.bincount1 is the preferrably used when different
scaleBins <- function(start, end, target.bincount1 = 5, target.bincount2 = 5) {
  if(end < start) {
    # switch
    foo <- start; start <- end; end <- foo
  }
  if(end - start >= max(target.bincount1, target.bincount2)) {
    # make some good divisions to produce 4 or 5 bins
    # the last one always works but may be not so well-suited
    bins.trial <- list(
        seq(start, end, length.out = target.bincount1), 
        seq(start, end, length.out = target.bincount2), 
        seq(start +1, end, length.out = target.bincount1), 
        seq(start +1, end, length.out = target.bincount2), 
        seq(start +2, end, length.out = target.bincount1), 
        seq(start +2, end, length.out = target.bincount2),
        seq(start +3, end, length.out = target.bincount1), 
        seq(start +3, end, length.out = target.bincount2), 
        round(seq(start +4, end, length.out = 5))
    )
    # which trials are whole-numbered? (the last one always will be -> termination)
    sel.idx <- which(sapply(bins.trial, function(bins) all(bins-round(bins) == 0)))
    # use the first of them
    return(unlist(bins.trial[min(sel.idx)]))
  } else {
    return(start:end)
  }
}

