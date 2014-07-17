manhattanplot <- function(
                    gwas.dataset, 
                    highlight.logp = c(6, 7.3), 
                    highlight.win = 125000, 
                    highlight.color = NULL, 
                    highlight.text = c("SNP", "genes"),
                    highlight.cex = c(0.7, 1),
                    highlight.fontface = c("italic", "bold"),
                    highlight.lines = FALSE,
                    ticks.y = FALSE, 
                    max.y = NULL,
                    reduce.dataset = TRUE, 
                    plot.title = NULL, 
                    biomart.config = biomartConfigs$hsapiens, 
                    use.buffer = FALSE, 
                    toFile = nextFilename("manhattanplot", "pdf")
                  ) {

  ######### prepare data #########

  if(is.null(gwas.dataset) || missing(gwas.dataset))
    stop("Argument gwas.dataset has to be provided.")
  
  gwasdat <- readGWASdatasets(gwas.dataset)
  
  if( !("SNP" %in% names(gwasdat) && "CHR" %in% names(gwasdat) && "BP" %in% names(gwasdat) && "P" %in% names(gwasdat) ) )
    stop("The GWAS dataset has to contain columns CHR, BP, and P.")

  if(is.null(highlight.logp) || length(highlight.logp) < 1)
    stop("Argument highlight.logp has to be set.\n")
  
  if(is.null(highlight.color)) {
    highlight.color <- rev(rainbow(length(highlight.logp)))
  } else if(!all(highlight.color %in% colors())) {
    stop("Highlighting colors not defined or unknown (see colors())")
  }
  
  if(is.null(highlight.cex)) {
    stop("Argument highlight.cex has to be set.\n")
  }

  if(is.null(highlight.fontface)) {
    stop("Argument highlight.fontface has to be set.\n")
  }
  
  if(is.null(highlight.text)) {
    highlight.text <- "none"
  }
  
  if(is.null(highlight.lines) || !is.logical(highlight.lines)) {
    stop("Argument highlight.lines has to be set and has to be logical.\n")
  }
  
  # preprocess data
  gwasdat <- vectorElements(gwasdat[, c("SNP", "CHR", "BP", "P")])
  gwasdat$P  <- as.numeric(gwasdat$P)
  gwasdat$BP <- as.numeric(gwasdat$BP)
  gwasdat <- na.omit(gwasdat)
  gwasdat$logp <- -log10(gwasdat$P)
  gwasdat <- gwasdat[gwasdat$P > 0 & gwasdat$P < 1, ]
  
  if(reduce.dataset > 0) {
#    gwasdat <- gwasdat[gwasdat$logp > runif(nrow(gwasdat))^2 | gwasdat$logp < 1e-04, ]
#    gwasdat <- gwasdat[gwasdat$logp > (runif(nrow(gwasdat))/1.1)^2 | gwasdat$logp < 1e-04, ]
    gwasdat <- reduceGWAS(gwasdat, reduce.dataset)
  }
  
  # check max.y argument (when not user defined, set to max(gwasdat$logp))
  if(is.null(max.y) || !is.numeric(max.y) || length(max.y) != 1 || max.y <= 0) {
    max.y <- max(gwasdat$logp)
  } else {
    # remove entries that exceed max.y
    gwasdat <- gwasdat[gwasdat$logp < max.y, ]
  }
  
  # translate base position to coords (-> gwasdat$pos), making chromosmes sequential
  
  # use mapping to generate numeric chromosomes so that we can sequentialize base positions
  chr.map  <- data.frame(CHR = chromosort(unique(gwasdat$CHR)), CHR.mapped = 1:length(unique(gwasdat$CHR)))
  gwasdat  <- merge(gwasdat, chr.map)
  chr.max.base <- tapply(gwasdat$BP, factor(gwasdat$CHR.mapped), max)  # for each chromosome the bp length (vector name = chrnames)
  chr.shift <- sapply(
                  names(chr.max.base), 
                  function(chr) sum(as.double(chr.max.base[as.numeric(names(chr.max.base)) <= as.numeric(chr)])) - chr.max.base[chr],
                  USE.NAMES = FALSE
                )
  for(i in names(chr.shift))
    gwasdat[gwasdat$CHR.mapped == i, "pos"] <- gwasdat[gwasdat$CHR.mapped == i, "BP"] + chr.shift[as.character(i)]

  # reorganize (and recycle) parameters when multiple p thresholds are given
  if(length(highlight.logp) <= 1) { # data frame stuff does not work with single row...
    highlight.logp <- highlight.logp[1] 
    highlight.win <- highlight.win[1]
    highlight.color <- highlight.color[1] 
    highlight.cex <- highlight.cex[1]
    highlight.fontface <- highlight.fontface[1] 
    highlight.text <- highlight.text[1]
    highlight.lines <- highlight.lines[1]
  } else {
    param <- cbind(highlight.logp, highlight.win, highlight.color, highlight.cex, highlight.fontface, highlight.text, highlight.lines)
    # in case that length(highlight.logp) < another argument vector, make logp thresholds unique
    param <- param[!duplicated(param[, "highlight.logp"]), ] 
    # param will always have > 1 rows (length(highlight.logp) > 1) and will therefore be still a matrix / data frame
    # highlight tresholds have to be sorted
    param <- param[order(highlight.logp), ]
    highlight.logp <- param[, "highlight.logp"]
    highlight.win <- param[, "highlight.win"]
    highlight.color <- param[, "highlight.color"]
    highlight.cex <- param[, "highlight.cex"]
    highlight.fontface <- param[, "highlight.fontface"]
    highlight.text <- param[, "highlight.text"]
    highlight.lines <- param[, "highlight.lines"]
  }
  
  ######### set highlighted SNP color and shape #########
  
  message("Identifying highlighted regions...")
  clr <- data.frame(CHR.mapped = sort(unique(as.vector(gwasdat$CHR.mapped))))
  clr[seq(1, nrow(clr), 2), "color"] <- "grey10"
  if(nrow(clr) > 1)
    clr[seq(2, nrow(clr), 2), "color"] <- "grey50"
  gwasdat <- merge(gwasdat, clr)
  
  gwasdat$x <- gwasdat$pos / max(gwasdat$pos)
  gwasdat$y <- gwasdat$logp / max.y
  highlight.lines.y <- as.numeric(highlight.logp) / max.y
  
  # highl only contains the PEAK data point of highlighted regions
  highl <- mapply(
    function(logp, win, col, cex, font, text) {
      highl    <- gwasdat[gwasdat$logp > as.numeric(logp), ]
      highl    <- removeNeighborSnps(highl, as.numeric(win))
      highl$CHR.mapped <- as.numeric(as.vector(highl$CHR.mapped))
      highl$BP <- as.numeric(as.vector(highl$BP))
      highl$P <- as.numeric(as.vector(highl$P))
      highl$color <- rep(col, nrow(highl))
      highl$cex <- rep(cex, nrow(highl))
      highl$font <- rep(font, nrow(highl))
      highl$text <- rep(text, nrow(highl))
      return(highl)
    }, 
    highlight.logp, 
    highlight.win, 
    highlight.color,
    highlight.cex,
    highlight.fontface,
    highlight.text, 
    SIMPLIFY = FALSE
  )
  # for multiple p threshs: make real intervals with lower bounds 
  # (currently only upper bound exists for each thresh, i.e. > logp)
  if(length(highlight.logp) > 1)
    for(i in 1:(length(highlight.logp)-1))
      highl[[i]] <- highl[[i]][as.numeric(highl[[i]]$logp) <= as.numeric(highlight.logp[i+1]), ]
  

  # shape of all ordinaty data points: circle
  gwasdat$shape <- 20
  
  # all highlighted data points: triangle
  # also set color and cex of highlighted data points
  mapply(
    function(highl, win, color, cex) {
      if(nrow(highl) > 0) {
        for(idx in 1:(nrow(highl))) {
          gwasdat.highl.idx <- which(
                                 as.numeric(gwasdat$CHR.mapped) == as.numeric(highl[idx, "CHR.mapped"]) & 
                                 gwasdat$BP >= as.numeric(highl[idx, "BP"]) - win & 
                                 gwasdat$BP <= as.numeric(highl[idx, "BP"]) + win
                               )
          gwasdat[gwasdat.highl.idx, "color"] <<- color
          gwasdat[gwasdat.highl.idx, "shape"] <<- 2
          gwasdat[gwasdat.highl.idx, "cex"] <<- cex
        }
      }
    },
    highl, 
    as.numeric(highlight.win), 
    highlight.color,
    highlight.cex
  )
  
  
  ######### plot points and genes #########
  if(!is.null(toFile) && is.character(toFile) && length(toFile) == 1)
    pdf(toFile, width = 4.8, height = 2, pointsize = 4)
  
    
  # use npc: everything scaled to papersize percent
  message("Scatterplot...")
  # main plot area size (leaves 0.125 space up and down and 0.08 left for axes viewports)
  vp.plot <- viewport(x = 0.52, width = 0.88, height = 0.75, name = "plot")
  pushViewport(vp.plot)
  
  gwasdat.plain <- gwasdat[gwasdat$shape == 20, ]
  if(nrow(gwasdat.plain) > 0) {
    grid.points(
        x = unit(gwasdat.plain$x, "npc"), 
        y = unit(gwasdat.plain$y, "npc"),
        pch = gwasdat.plain$shape,
        size = unit(0.0004, "npc"),
        gp = gpar(col = gwasdat.plain$color)
    )
  } else {
    warning("Manhattanplot: Source file does only contain highlighted data")
  }
  
  # plot calibration lines if desired
  highlight.lines <- as.logical(highlight.lines) 
  if(any(highlight.lines)) {
    grid.polyline(
        x = unit(rep(0:1, sum(highlight.lines)), "npc"), 
        y = unit(rep(highlight.lines.y, each = sum(highlight.lines)), "npc"),
        id = rep(which(highlight.lines), each = 2), 
        gp = gpar(col = highlight.color, lwd = 0.8 * as.numeric(highlight.cex))
    )
  }

  gwasdat.highl <- gwasdat[gwasdat$shape != 20, ]
  if(nrow(gwasdat.highl) > 0) {
    
    grid.points(
      x = unit(gwasdat.highl$x, "npc"), 
      y = unit(gwasdat.highl$y, "npc"),
      pch = 25,
      size = unit(0.004 * as.numeric(gwasdat.highl$cex)^1.7, "npc"),
      gp = gpar(col = gwasdat.highl$color, fill = gwasdat.highl$color)
    )

    # text (gene) annot
    highl.df <- list2df(highl) 
    if(nrow(highl.df) > 0) {

      if(any(highl.df$text == "SNP")) {
        highl.snp <- highl.df[highl.df$text == "SNP", ]
        grid.text(
            label = highl.snp$SNP, 
            x = unit(highl.snp$x, "npc"), 
            y = unit(highl.snp$y, "npc") + unit(0.013 * as.numeric(highl.snp$cex), "npc"),
            just = "bottom",
            gp = gpar(col = highl.snp$color, fill = highl.snp$color, cex = 0.7 * as.numeric(highl.snp$cex), fontface = highl.snp$font)
        )
      }
      
      if(any(highl.df$text == "genes")) {
        message("Geneplot...")
        highl.gene.toannot <- highl.df[highl.df$text == "genes", ]
        colnames(highl.gene.toannot)[colnames(highl.gene.toannot) == "CHR.mapped"] <- "chrmp"
        highl.gene <- snp2gene.prox(snps = highl.gene.toannot, by.genename = TRUE, level = 1, biomart.config = biomart.config, use.buffer = use.buffer)
        
        # check SNPs without gene annotation, plot SNP id on top then 
        highl.gene.noannot <- highl.gene.toannot[!(highl.gene.toannot$SNP %in% unique(highl.gene$SNP)), ]
        if(nrow(highl.gene.noannot) > 0)
          grid.text(
              label = highl.gene.noannot$SNP, 
              x = unit(highl.gene.noannot$x, "npc"), 
              y = unit(highl.gene.noannot$y, "npc") + unit(0.0175 * as.numeric(highl.gene.noannot$cex), "npc"),
              just = "bottom",
              gp = gpar(col = highl.gene.noannot$color, fill = highl.gene.noannot$color, cex = 0.7 * as.numeric(highl.gene.noannot$cex), fontface = highl.gene.noannot$font)
          )
        
        if(nrow(highl.gene) > 0) {
          
          highl.gene.left <- highl.gene[highl.gene$direction == "up", ]
          if(nrow(highl.gene.left) > 0)
            grid.text(
                label = highl.gene.left$genename, 
                x = unit(highl.gene.left$x, "npc") - unit(0.005 * as.numeric(highl.gene.left$cex)^2, "npc"), 
                y = unit(highl.gene.left$y, "npc"),
                just = "right",
                gp = gpar(col = highl.gene.left$color, fill = highl.gene.left$color, cex = 0.7 * as.numeric(highl.gene.left$cex), fontface = highl.gene.left$font)
            )
          
          highl.gene.right <- highl.gene[highl.gene$direction == "down", ]
          if(nrow(highl.gene.right) > 0)
            grid.text(
                label = highl.gene.right$genename, 
                x = unit(highl.gene.right$x, "npc") + unit(0.005 * as.numeric(highl.gene.right$cex)^2, "npc"), 
                y = unit(highl.gene.right$y, "npc"),
                just = "left", 
                gp = gpar(col = highl.gene.right$color, fill = highl.gene.right$color, cex = 0.7 * as.numeric(highl.gene.right$cex), fontface = highl.gene.right$font)
            )
          
          highl.gene.top <- highl.gene[highl.gene$direction == "cover", ]
          if(nrow(highl.gene.top) > 0)
            grid.text(
                label = highl.gene.top$genename, 
                x = unit(highl.gene.top$x, "npc"), 
                y = unit(highl.gene.top$y, "npc") + unit(0.0175 * as.numeric(highl.gene.top$cex), "npc"),
                just = "bottom",
                gp = gpar(col = highl.gene.top$color, fill = highl.gene.top$color, cex = 0.7 * as.numeric(highl.gene.top$cex), fontface = highl.gene.top$font)
            )
        }
      }
    }
  }
  
  ######### plot axes #########

  message("Adding axes...")
  # y axis
  vp.yaxis <- viewport(x = unit(0.04, "npc"), width = 0.08, height = 0.75)
  upViewport()
  pushViewport(vp.yaxis)
  
  y.maxtick <- floor(max.y)
  # at most 8 ticks
  ticks <- seq(1, y.maxtick, by = ceiling(y.maxtick/8))
  # add y axis numbers
  grid.text(
    label = ticks, 
    x = unit(0.6, "npc"), 
    y = ticks / max.y,
    gp = gpar(col = "grey30", cex = 0.8)
  )
  if(ticks.y) {
    # add y axis line with ticks
    grid.polyline(
        x = rep(c(0.75, 0.85), each = length(ticks)), 
        y = rep(ticks / max.y, 2),
        id = rep(1:length(ticks), 2),
        gp = gpar(lwd = 0.7)
    )
    grid.lines(
        x = c(0.85, 0.85), 
        y = c(min(ticks / max.y), max(ticks / max.y)), 
        gp = gpar(lwd = 0.7)
    )
  }
  
  grid.text(
    label = "-log10(p)", 
    x = unit(0.3, "npc"), 
    y = unit(0.5, "npc"),
    rot = 90, 
    gp = gpar(col = "grey30", fontface = "bold", cex = 0.8)
  )
  
  # x axis
  vp.axes <- viewport(x = 0.52, y = unit(0.05, "npc"), width = 0.88, height = 0.1)
  upViewport()
  pushViewport(vp.axes)
  
  ticks.pos <- NULL
  for(i in chr.map$CHR.mapped) 
    ticks.pos = c(ticks.pos, mean(gwasdat[gwasdat$CHR.mapped == i, "x"]))
  
  grid.text(
    label = chr.map$CHR[odd(1:nrow(chr.map))], 
    x = ticks.pos[odd(1:nrow(chr.map))], 
    y = unit(0.66, "npc"),
    gp = gpar(col = "grey30", cex = 0.7)
  )
  
  if(nrow(chr.map) > 1)
    grid.text(
      label = chr.map$CHR[even(1:nrow(chr.map))], 
      x = ticks.pos[even(1:nrow(chr.map))], 
      y = unit(0.775, "npc"),
      gp = gpar(col = "grey30", cex = 0.7)
    )
  
  grid.text(
    label = "Chromosome", 
    x = unit(0.5, "npc"), 
    y = unit(0.28, "npc"),
    gp = gpar(col = "grey30", fontface = "bold", cex = 0.8)
  )
  
  
  ######### plot title #########
  
  if(!is.null(plot.title)) {
    vp.title <- viewport(y = unit(0.975, "npc"), width = 0.9, height = 0.05)
    upViewport()
    pushViewport(vp.title)
    
    grid.text(
      label = plot.title, 
      x = unit(0.5, "npc"), 
      y = unit(0.5, "npc"),
      gp = gpar(col = "grey30", cex = 0.8)
    )
  }

  if(!is.null(toFile) && is.character(toFile) && length(toFile) == 1)
    dev.off()
  message("Done.")
  return(NULL)
}


# For chromosome identifiers, numbers and characters are often mixed. 
# This function takes a vector of numeric or character values and 
# returns a properly ordered character vector.  
chromosort <- function(chromnames) {
  only.num  <- suppressWarnings(na.omit(as.numeric(chromnames)))
  only.char <- suppressWarnings(chromnames[is.na(as.numeric(chromnames))])
  return(c(
          only.num[order(only.num)], 
          only.char[order(only.char)]
  ))
}


# gwasdat needs to be a data frame with a column P of p-values
# level is an integer, the higher, the smaller the returned values. May not be negative.
# returns the pruned gwasdat data frame
reduceGWAS <- function(gwasdat, level) {
  
  # helper function 
  # generates a vector of random numbers 
  # this can be used to prune the argument vector pvals
  prunePs <- function(pvals, level, min = quantile(pvals, 0.002), max = quantile(pvals, 0.6)) {
    if(level <= 0) {
      rep(1, length(pvals))
    } else {
      pmin(
          runif(length(pvals), min = min, max = max), 
          prunePs(pvals, level -1, min, max)
      )
    }
  }
  
  remove.idx <- gwasdat$P > prunePs(gwasdat$P, level)
  gwasdat[!remove.idx, ] 
}
