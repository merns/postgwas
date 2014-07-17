gwas2network.plot <- function(
                   g = g,
                   layout = layout.fruchterman.reingold(g, area = vcount(g)^3),
                   legend.autolayout = TRUE, 
                   vertex.to.edge.space = 1.3, 
                   device = pdf, 
                   ...
                 ) {

  # general hint: The vertex IDs serve as indices of rows in the layout matrix (i.e. vertex coordinates)
  # vertex IDs can be obtained by stating as.vector(V(graph))
  
  if(!(class(g) == "igraph") || vcount(g) < 1)
    stop("Argument g is not an object of class 'igraph' or empty.")
                   
	# there is a legacy package 'igraph0' that masks igraph and leads to errors
	# ensure that our igraph is first on search path 
	# and remember its current position for restore at function exit
#	igraph.pos.orig <- which(search() == "package:igraph")
#	suppressWarnings(detach("package:igraph", force = TRUE))
#	suppressWarnings(suppressPackageStartupMessages(library(igraph, pos = 2, quietly = TRUE)))

  args <- list(...)
  if(any(c("width", "height", "units") %in% names(args)))
    stop("Arguments 'width', 'height' and 'units' cannot be specified in this function\n")
  if(!is.null(args$v.weight.max)) {
    vwm <- args$v.weight.max
    args$v.weight.max <- NULL
  } else {
    vwm <- max(V(g)$weight)    
  }
  pointfac <- ifelse(is.null(args$pointsize), 12, args$pointsize) / 12
  V(g)$radius <- V(g)$weight / vwm
  # set vertex frames invisible (last two characters 00 = transparent, FF solid)
  V(g)$frame.color <- "#00000000"
  V(g)$label.fontface <- "1" # plain
  
  # use x/y coordinates starting from 0 
  layout[, 1] <- layout[, 1] - min(layout[, 1])
  layout[, 2] <- layout[, 2] - min(layout[, 2])
  
  # sizes in inches
  mainvp.size <- (max(layout))^0.45
  mainvp.size <- if(mainvp.size < 10) 10 else mainvp.size
  mainvp.size <- mainvp.size * (if(!is.null(E(g)$label)) 1.25 else 1.1)  # add some (border) space for (large) labels
  scalevp.size <- list(x = -1.5, y = 0)
  enrichvp.size <- list(x = -1.5, y = 0)
  if(legend.autolayout) {
    # reserve two viewports for scale and enrichment labels at the bottom left and right
    # will plot scale and legend vertices one above each other (by vertex weight) in these viewports
    # find indices of these vertices
    scale.idx <- grepl("postgwasScale", V(g)$name, fixed = TRUE)
    enrich.idx <- grepl("postgwasEnrich", V(g)$name, fixed = TRUE)
    # introduce linebreaks in long labes (e.g. GO-term descriptions)
    if(any(enrich.idx))
      V(g)$label[enrich.idx] <- V(g)$label[enrich.idx]
    # we have to open a dummy device to estimate the width of strings
    if(identical(device, png) || identical(device, jpeg) || identical(device, tiff) || identical(device, bmp)) 
      do.call(device, c(width = 1, height = 1, units = "in", args))
    else
      do.call(device, c(width = 1, height = 1, args))

    # determine required viewport sizes
    if(any(enrich.idx)) {
      enrich.labelw <- max(strwidth(V(g)$label[enrich.idx], units = "inches"))
      enrich.labelh <- rev(strheight(V(g)$label[enrich.idx], units = "inches"))
      enrichvp.size$x <- enrich.labelw + 0.8 * pointfac
      enrichvp.size$y <- sum(enrich.labelh) + (sum(enrich.idx) +1) * 0.4 * pointfac
      enrich.ypos <- rev(sapply(
          1:length(enrich.labelh), 
          function(idx) sum(enrich.labelh[1:idx] + 0.4 * pointfac)
      ))
    }
    
    if(any(scale.idx)) {
      scale.labelw <- max(strwidth(V(g)$label[scale.idx], units = "inches"))
      scale.labelh <- max(strheight(V(g)$label[scale.idx], units = "inches"))
      # vertices sizes are normalized and can be assumed constant (except labelw, and pointsize alterations)
      scalevp.size$x <- 1.25 * pointfac + scale.labelw 
      scalevp.size$y <- 0.75 * sum(scale.idx) + scale.labelh
    }
    dev.off() # close the dummy device
  }
  
  layout.normed <- layout / max(layout)
  layout.normed[is.na(layout.normed)] <- 0
  layout.normed <- as.data.frame(layout.normed) # we need that because subsetting a matrix to one row will yield a vector!

  # the size of the output device is 
  # main + legend + 1.5 inches outer border at both sides (for long labels) 
  # and two times 1.5 inch space on the x axis between main and legend viewports
  papersize <- list(
      x = 1.5 + scalevp.size$x + 1.5 + mainvp.size + 1.5 + enrichvp.size$x +1.5, 
      y = 1.5 + max(mainvp.size, enrichvp.size$y, scalevp.size$y) + 1.5
  )
  if(identical(device, png) || identical(device, jpeg) || identical(device, tiff) || identical(device, bmp)) 
    do.call(device, c(width = papersize$x, height = papersize$y, units = "in", args))
  else
    do.call(device, c(width = papersize$x, height = papersize$y, args))
  
  # draw legend separately, if applicable
  if(legend.autolayout) {
    
    if(any(enrich.idx)) {
      pushViewport(viewport(
              x = unit(1, "npc") - unit(1.5, "inches"), y = 1.5, 
              width = enrichvp.size$x, 
              height = enrichvp.size$y, 
              default.units = "inches", 
              just = c("right", "bottom")
          ))
      grid.rect()
      grid.text(
          label = V(g)[enrich.idx]$label, 
          x = 0.4 / enrichvp.size$x, # in npc, 0.4 inches to both sides (labels are left-justified) 
          y = enrich.ypos / enrichvp.size$y,
          just = c("left", "top"),
          gp = gpar(col = V(g)[enrich.idx]$label.color)
      )
      upViewport()
    }
    
    if(any(scale.idx)) {
      pushViewport(viewport(
              x = 1.5, y = 1.5,  
              width = scalevp.size$x, 
              height = scalevp.size$y,
              default.units = "inches", 
              just = c("left", "bottom")
          ))
      grid.rect()
      # equally distributed, in npc units
      scale.ypos <- seq(0, 1, length.out = sum(scale.idx)+2)[2:(sum(scale.idx)+1)]  - scale.labelh / scalevp.size$y
      scale.xpos <- 0.5 - scale.labelw / scalevp.size$y
      layout.normed[as.vector(V(g)[scale.idx]), ] <- as.matrix(data.frame(scale.xpos, scale.ypos))
      plotVertexshapes(V(g)[scale.idx], layout = layout.normed)
      plotVertexlabels(V(g)[scale.idx], layout = layout.normed)
      upViewport()
    }
    
    # remove the already plotted legend vertices from the main graph
    g <- g - V(g)[scale.idx | enrich.idx]
  }
  

  # draw main plot (x and y centered)
  main.vp <- viewport(
      x = 1.5 + scalevp.size$x + 1.5 + mainvp.size / 2, 
      width = mainvp.size, 
      height = mainvp.size, 
      default.units = "inches"
  )
  pushViewport(main.vp)

  # plot edges first (are in background)
  if(ecount(g) > 0) {
    el <- na.omit(get.edgelist(g, names = F))   # vertex indices, eg V(g)$name[el], and layout index
    lapply(
      1:nrow(el), 
      function(idx) {
        # we define a start and end vertex for the edge
        x <- c(start = layout.normed[el[idx, 1], 1], end = layout.normed[el[idx, 2], 1])
        y <- c(start = layout.normed[el[idx, 1], 2], end = layout.normed[el[idx, 2], 2])
        
        # edges should be drawn from the vertex circle border
        # determine the point on the circle border where the edge should start
        # this is the intersection of edge and circle, where the edge is the hypothenusis or a triangle
    
        # radius of both circles (we add some space so that edge does not start directly at vertex border)
        r.start <- V(g)$radius[el[idx, 1]] * vertex.to.edge.space
        r.stop  <- V(g)$radius[el[idx, 2]] * vertex.to.edge.space
    
        # triangle
        adjacent <- abs(x["end"] - x["start"])
        opposite <- abs(y["end"] - y["start"])
        alpha    <- atan(opposite / adjacent)
        # for the triangle configuration with start vertex to the left and below end vertex
        x.adj <- cos(alpha)
        y.adj <- sin(alpha)

        if(x["start"] > x["end"]) 
          x.adj <- -x.adj
        if(y["start"] > y["end"])
          y.adj <- -y.adj
        
        x.just <- c(start = x.adj * r.start, end = -x.adj * r.stop)
        y.just <- c(start = y.adj * r.start, end = -y.adj * r.stop)
        grid.lines(
          x  = unit(x, "npc") + unit(x.just, "char"), 
          y  = unit(y, "npc") + unit(y.just, "char"),
          gp = gpar(
                 col = E(g)$color[idx]
               )
        )
        
        # draw labels
        if(!is.null(E(g)$label)) {
          grid.text(
            splitSentence(E(g)$label[idx]), 
            x  = unit((x["start"] + x["end"]) / 2, "npc"), 
            y  = unit((y["start"] + y["end"]) / 2, "npc"),
            gp = gpar(
                   col = E(g)$color[idx], 
                   fontface = "italic", 
                   cex = 0.7  # make edge labels a bit smaller
                  )
          )
        }
      }
    )
  }
  
  # plot vertices
  if(vcount(g) > 0) {

    # category information is optional: all.cat can be empty, then category will be set to 'circle' globally
    # determine all categories used. We have to see that from column names, because
    # when for one category all vertices are a complete subset of another categories vertices, 
    # this will not be listed in category.main.
    all.cat <- list.vertex.attributes(g)[list.vertex.attributes(g) %in% c("circle", "square", "diamond")]
    
    if(length(all.cat) > 3)
      stop("Too many categories - cannot plot that.\n")
    
    if(length(all.cat) <= 0) {
      V(g)$category.main <- "circle"
    } else {
      # set unknown category values
      V(g)$category.main[is.na(V(g)$category.main)] <- all.cat[1]
      V(g)$category.main[!(V(g)$category.main %in% all.cat)] <- all.cat[1]
      for(curr.cat in all.cat)
        set.vertex.attribute(g, curr.cat, value = na.set(get.vertex.attribute(g, curr.cat)))
    }
    
    # plot all vertices with the shape of their main category
    for(curr.cat in all.cat) {
      # sequence of vertices for the current category to plot
      vs.curr.cat <- V(g)[V(g)$category.main == curr.cat]
      plotVertexshapes(vs.curr.cat, layout = layout.normed, type = curr.cat)
    }
    
    # now plot multi-category vertices on top of the main category vertices
    for(curr.cat in all.cat) {
      # get sequence of vertices having the current category as alternative category
      # by na.set (to FALSE), ignores vertices with NA multi category setting
      multi.idx <- na.set(get.vertex.attribute(g, curr.cat))
      g2 <- g
      # we plot only the border of the shape of the alternative category (last two characters 00 = transparent, FF solid)
      V(g2)[multi.idx]$color <- "#00000000"
      # V(subg)$color <- adjustcolor(vs.multicat$color, red.f = 0.9, blue.f = 0.9, green.f = 0.9, offset = c(0, 0, 0, 1))
      V(g2)[multi.idx]$frame.color <- "white"
      V(g2)[multi.idx]$radius <- get.vertex.attribute(g, paste(curr.cat, "weight", sep = "."))[multi.idx] / vwm
      V(g2)[multi.idx]$marked <- get.vertex.attribute(g, paste(curr.cat, "marked", sep = "."))[multi.idx]
      plotVertexshapes(V(g2)[multi.idx], layout = layout.normed, type = curr.cat)
      rm(g2)
    }

    # plot vertex text (is always drawn once per gene, at the main category vertex)
    # set different font for multi category
    for(curr.cat in all.cat)
      V(g)$label.fontface[get.vertex.attribute(g, curr.cat)] <- "4" # bold-italic
    
    plotVertexlabels(V(g), layout.normed)

  }
  
  
  if(names(dev.cur()) %in% c('X11', 'windows', 'win.graph')) {
    # on-screen device was used, keep open
  } else {
    dev.off()
  }
    
  # restore search path
#  suppressWarnings(detach("package:igraph", force = TRUE))
#  suppressWarnings(suppressPackageStartupMessages(library(igraph, pos = igraph.pos.orig, quietly = TRUE)))

}



# v is an igraph vertex sequence (e.g. returned by V(graph)) having attributes radius, color, frame.color, marked
# the vertex IDs have to be indices for rows of the layout matrix
# layout may contain additional coordinates of other vertices (more rows)
# layout is in npc
plotVertexshapes <- function(v, layout, type = "circle") {
  
  if(!(class(v) == "igraph.vs") || is.null(v$radius) || is.null(v$color) || is.null(v$frame.color) || is.null(v$marked))
    stop("plotVertexshapes: Argument v is not a vertex sequence or lacks one of the attributes radius, color, frame.color, marked.\n")
  
  if(!is.null(v) && length(v) > 0) {
    
    if(max(v) > nrow(layout) || min(v) < 1) # max and min select max/min vertex indices
      stop("plotVertexshapes: Vertex IDs do not index the layout coordinate matrix rows.\n")
    
    # only layout rows of vertex sequence
    layout <- layout[as.vector(v), ]
    
    if(type == "circle") {
      grid.circle(
        x = layout[, 1], 
        y = layout[, 2], 
        r = unit(v$radius, "char"),
        gp = gpar(
          lex = 4 * v$radius,
          col = as.vector(v$frame.color),
          fill = as.vector(v$color)
        )
      )
    }

    if(type == "square") {
      grid.rect(
        x = layout[, 1], 
        y = layout[, 2], 
        width = unit(sqrt(pi * v$radius^2), "char"), 
        height = unit(sqrt(pi * v$radius^2), "char"), 
        gp = gpar(
          lex = 4 * v$radius,
          col = as.vector(v$frame.color),
          fill = as.vector(v$color)
        )
      )
    }

    if(type == "diamond") {
      grid.polygon(
        x = unit.c(   # x coords of the four corners of the diamond
              unit(layout[, 1], "npc"), 
              unit(layout[, 1], "npc") + unit(sqrt(pi * v$radius^2) /2, "char"), 
              unit(layout[, 1], "npc"), 
              unit(layout[, 1], "npc") - unit(sqrt(pi * v$radius^2) /2, "char")
            ), 
        y = unit.c(
              unit(layout[, 2], "npc") + unit(sqrt(pi * v$radius^2) /2, "char"), 
              unit(layout[, 2], "npc"), 
              unit(layout[, 2], "npc") - unit(sqrt(pi * v$radius^2) /2, "char"), 
              unit(layout[, 2], "npc")
            ),
        id = rep(1:nrow(v), 4), 
        gp = gpar(
          lex = 4 * v$radius,
          col = as.vector(v$frame.color),
          fill = as.vector(v$color)
          )
        )
    }

    # plot dual annot marks on vertex positions
    if(!is.null(v$marked) && sum(v$marked) > 0) {
      grid.points(
        x = layout[v$marked, 1], 
        y = layout[v$marked, 2], 
        size = unit(v$radius[v$marked], "char"),
        pch = 4,
        gp = gpar(
          lex = 2.5 * v$radius[v$marked], 
          col = "white"
        )
      )
    }
  }
  
}


# v is an igraph vertex sequence (e.g. returned by V(graph)) having attributes label, radius, label.color, label.fontface
# the vertex IDs have to be indices for rows of the layout matrix
# layout may contain additional coordinates of other vertices (more rows)
# layout is in npc
plotVertexlabels <- function(v, layout) {
  
  if(!(class(v) == "igraph.vs") || is.null(v$label) || is.null(v$radius) || is.null(v$label.color) || is.null(v$label.fontface))
    stop("plotVertexlabels: Argument v is not a vertex sequence or lacks one of the attributes label, radius, label.color, label.fontface.\n")
  
    if(!is.null(v) && length(v) > 0) {
      
      if(max(v) > nrow(layout) || min(v) < 1) # max and min select max/min vertex indices
        stop("plotVertexshapes: Vertex IDs do not index the layout coordinate matrix rows.\n")
    
      layout <- layout[as.vector(v), ]
      
      grid.text(
        label = v$label, 
        x = unit(layout[, 1], "npc") + unit(v$radius, "char"), 
        y = unit(layout[, 2], "npc") + unit(v$radius, "char"),
        just = "left",
        gp = gpar(
          col = v$label.color, 
          fontface = as.numeric(as.vector(v$label.fontface))
        )
      )
  }
}
