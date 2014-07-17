regionalplot <- function(
    snps, 
    gwas.datasets, 
    window.size = 1000000, 
    biomart.config = biomartConfigs$hsapiens,
    use.buffer = FALSE,
    plot.genes = TRUE, 
    draw.snpname = "auto",
    ld.options = list(
        gts.source = 2, 
        max.snps.per.window = 200, 
        rsquare.min = 0.2, 
        show.rsquare.text = FALSE),
    var.options = NULL,
    out.format = list(
        file = "pdf", 
        panels.per.page = "auto",
        global.scale = "auto", 
        paper.height = 11.7, 
        paper.width = 8.3),
    max.logp = FALSE, 
    cores = 1,
    ytracks = ytracks.regionalplot, 
    panel.add = function(...) return(NULL)
) {

  # avoid unknown binding warnings in codecheck tool
  pheno <- "NULL"
  
  if(missing(gwas.datasets) || is.null(gwas.datasets)) 
   stop("Function parameter 'gwas.datasets' needs to be specified")
 
  if(is.null(biomart.config) || !is.list(biomart.config))
    stop("Argument 'biomart.config' has to be a list")
  
  if(missing(snps) || is.null(snps)) 
   stop("Function parameter 'snps' needs to be specified")

  if(!is.data.frame(snps) || nrow(snps) < 1) 
    stop("Function parameter snps is not a data frame or empty")

  if(is.null(snps$SNP)) 
    stop("Function parameter 'snps' does not contain 'SNP' column")
  
  if(is.null(ld.options)) {
    plot.ld <- FALSE
  } else {
    plot.ld <- TRUE
    # when missing, use default setting
    if(is.null(ld.options$rsquare.min)) ld.options$rsquare.min <- 0.2
    if(is.null(ld.options$max.snps.per.window)) ld.options$max.snps.per.window <- 200
    if(is.null(ld.options$show.rsquare.text)) ld.options$show.rsquare.text <- FALSE
    if(is.null(ld.options$gts.source)) 
      stop("Function parameter ld.options lacks an elements 'gts.source' defining the source for genotype retrieval.")
  }
  
  if(is.null(var.options)) {
    plot.variants <- FALSE
  } else {
    plot.variants <- TRUE
    if(!'VariantAnnotation' %in% installed.packages()[, 'Package']) {
      if(interactive() && readline("Package 'VariantAnnotation ' is not installed. Try to install? [Y/N] ") %in% c("Y", "y")) {
        biocLite <- NULL
        source("http://bioconductor.org/biocLite.R")
        biocLite("VariantAnnotation")
      } else {
        stop("Package VariantAnnotation is not installed but required when var.options are set.\n")		
      }
    }
    suppressPackageStartupMessages(stopifnot(require(VariantAnnotation, quietly = TRUE)))
    if(is.null(var.options$details)) var.options$details <- 3
    if(is.null(var.options$vcf.af.prune)) var.options$vcf.af.prune <- c(AF=50)
    if(length(var.options$vcf.af.prune) > 1) 
      stop("Function parameter var.options$vcf.af.prune has more that one element\n")
    if(is.null(var.options$hist.alpha)) var.options$hist.alpha <- 0.4
     if(is.null(var.options$vcf) || !is.list(var.options$vcf)) 
       stop("Function parameter var.options either misses 'vcf' component or it is not a list")
    # replace missing components with defaults
    if(!is.null(var.options$vcf$file))
      var.options$vcf <- list(var.options$vcf)
    var.options$vcf <- lapply(
      var.options$vcf,
      function(vcflist) {
        return(list(
          file            = if(is.null(vcflist$file) || !is.character(vcflist$file) || length(vcflist$file) > 1)
                              stop("Function parameter var.options misses 'vcf$file' component or has the wrong format.\n")
                            else
                              vcflist$file,
          remap.positions = if(is.null(vcflist$remap.positions)) 
                              stop("Function parameter var.options misses 'remap.positions' component")
                            else
                              vcflist$remap.positions, 
          chrom.map       = if(is.null(vcflist$chrom.map) || (is.data.frame(vcflist$chrom.map) && all(colnames(vcflist$chrom.map) %in% c("CHR", "CHR.VCF"))))
                              vcflist$chrom.map
                            else
                              stop("Function parameter var.options$vcf$chrom.map has the wrong format")
        ))
      }
    )
    if(length(var.options$vcf) == 1)
      var.options$vcf <- var.options$vcf[[1]]
  }

  if(is.null(out.format)) out.format <- list(file = NULL)
  if(is.null(out.format$paper.height)) out.format$paper.height <- 11.7
  if(is.null(out.format$paper.width)) out.format$paper.width <- 8.3
  if(is.null(out.format$panels.per.page) || out.format$panels.per.page == "auto") {
      out.format$panels.per.page <- 7 - sum(plot.genes, plot.ld, plot.variants)
  }
  if(is.null(out.format$global.scale) || out.format$global.scale == "auto") {
    # scaling factor is 1 when paper is A4 portrait
    global.scale <- sqrt(out.format$paper.height * out.format$paper.width / 97)
  } else {
    global.scale <- out.format$global.scale
  }

  if(cores > 1) {
	  if(!'parallel' %in% installed.packages()[, 'Package'])
		  stop("Function parameter cores > 1 requires the package 'parallel' to work correctly.\n")
	  suppressPackageStartupMessages(stopifnot(require(parallel, quietly = TRUE)))
  }
  
  # lazy initialization - assigned on use (thus buffer data can be used without connection)
  delayedAssign("snpds", bm.init.snps(biomart.config))
  delayedAssign("geneds", bm.init.genes(biomart.config))

  
  ##################### READ PVAL DATA AND EXTRACT REGIONS AROUND SNPS #####################
  
  p.all <- readGWASdatasets(gwas.datasets, assert.positions.match = TRUE)
  
  snps  <- data.regionalplot.snps(snps, p.all)
  if(nrow(snps) > 1000)
    stop("More than 1000 SNPs selected for plotting - cannot handle that.\n")

  # each region has at least one element - the centered SNP itself was checked to occur in pval file
  regions.df <- data.regionalplot.regions(snps, p.all, window.size)
  
  rm(p.all)
  
 
  ##################### BIOMART POSITIONS #####################

  # add BP.mapped column with positions from biomart (genes plotted have to match)
  if(plot.genes) {
    
    snps.bm <- bm.snps(config = biomart.config, ds = snpds, filter.val = regions.df$SNP, use.buffer= use.buffer)
    colnames(snps.bm)[colnames(snps.bm) == "chrom_start"] <- "BP.mapped"
    snps <- merge(snps, snps.bm, by.x = c("SNP", "CHR"), by.y = c("refsnp_id", "chrname"), all.x = TRUE)
    regions.df <- merge(regions.df, snps.bm, by.x = c("SNP", "CHR"), by.y = c("refsnp_id", "chrname"), all.x = TRUE)
    # impute missing positions chromosome-wise (if any)

    regions.df <- suppressWarnings(list2df(by(regions.df, regions.df$CHR, imputeRemappedPositions, remove.unmappable = TRUE)))
    snps <-  list2df(by(snps, snps$CHR, imputeRemappedPositions, remove.unmappable = TRUE))
    
    if(nrow(snps) < 1)
      stop("SNPs could not be mapped to biomart positions (were valid rs-numbers used? Do the chromosome names match biomart standards?)\n")
    
  } else {
    snps$BP.mapped <- snps$BP
    regions.df$BP.mapped <- regions.df$BP
    # missing positions do not occur, NAs have been removed from p.all
  }

  # add region start and end
  snps$region.startbp <- snps$BP.mapped - window.size/2
  snps$region.endbp <- snps$BP.mapped + window.size/2
  regions.df <- merge(regions.df, data.frame(packet.snp = snps$SNP, region.startbp = snps$region.startbp, region.endbp = snps$region.endbp))

  
  ##################### GET VARIANTS #####################
  
  if(plot.variants) {
    var.regions <- data.regionalplot.variants(snps, regions.df, var.options, biomart.config, snpds)
  }
    
  ##################### GET GENES #####################
  
  if(plot.genes) {
    genes.regions <- data.regionalplot.genes(
      snps, 
      regions.df, 
      window.size, 
      # half an inch of space for each gene label, measured in bp
      window.size / (out.format$paper.width * 2) * global.scale, 
      geneds, 
      biomart.config, 
      use.buffer)
  }
  
  ##################### GET GENOTYPES #####################

  if(plot.ld) {
    if(use.buffer && !is.null(get("ld.regionalplot", envir = postgwasBuffer))) {
      message("Using buffer data (ld)")
      ld.regions <- get("ld.regionalplot", envir = postgwasBuffer)
      if(!is.list(ld.regions) || !all(names(ld.regions) %in% snps$SNP))
        stop("Error in ld buffer data - clear buffer (run clearPostgwasBuffer()) or set use.buffer = FALSE")
    } else {
      ld.regions <- data.regionalplot.ld(regions.df, ld.options, cores)
      if(use.buffer) {
        assign("ld.regionalplot", ld.regions, envir = postgwasBuffer)
      }
    }
  }

  ##################### CALCUALTE SCALING, SPACING, LABELS FOR PLOT #####################
    
  # adjust default x-axis padding (additional space between graph and panel border)
  lattice.options.restore <- lattice.options()
  lattice.options.new <- list(axis.padding = list(numeric = 0))
  lattice.options(lattice.options.new)
    
  # calc y axis sections (space for pvalues, genes, ld)
  # ylim.upper: the largest value of data on the y-axis (some extra space will be added in the plot function)
  # ylim.space: a relative amount of space that can be used as universal unit for separating and borders
  # ylim.lower: the minimum value of the y-axis, this is artificial and depends on the following features drawn
  # genes.ysize: the amount of space on the negative y scale that is occupied by genes. Will vary by gene density
  # ld.ysize: the amount of space on the negative y scale that is occupied by LD
  # var.ysize: the amount of space on the negative y scale that is occupied by variants
    
  # p-values scale goes from 0 to ylim.upper 
  ylim.upper <- -1 * floor(log10(min(regions.df$P[!is.na(regions.df$P)])))
  if(max.logp != FALSE & ylim.upper > max.logp)
    ylim.upper <- max.logp
    
  # default separating and panel border space for y
  ylim.space <- ylim.upper / 100

  y.tracks <- ytracks(
                ylim.upper = ylim.upper, 
                ylim.space = ylim.space, 
                plot.genes = plot.genes, 
                plot.ld = plot.ld, 
                plot.variants = plot.variants, 
                var.options = var.options
              )
  ylim.lower <- - sum(y.tracks$yspace, y.tracks$ysize)

  # y axis labels
  # when ylim.upper is odd, add 1 (so there is a 'middle' position in 0:ylim.upper.odd)
  # y < 0: ticks based on number of tracks
  ylim.even <- ceiling(ylim.upper/2) *2
  y.at <- 0:ylim.even
  y.labels <- rep("   ", ylim.even +1)           # show labels for at most 5 ticks
  y.labels[seq(2, ylim.even, by = ceiling(ylim.even/5))] <- seq(1, ylim.even-1, by = ceiling(ylim.even/5))
  y.labels[ylim.even/2+1] <- paste( "-log10(p) ", y.labels[ylim.even/2+1])      # add "logp" at the median tick

  y.labels <- c(paste(y.tracks$name, collapse = "\n\n"), y.labels)
  y.at <- c(ylim.lower / 2, y.at)
  
  # page title & legend
  plot.title <- paste(
    "Regional association plot for multiple datasets", 
    if(plot.genes) "\nusing biomart positions", 
    if(plot.ld) 
      paste("\nLD options: ", 
            "genotypes =", if(class(ld.options$gts.source) == "snp.data") "custom"
               else if(is.numeric(ld.options$gts.source)) paste("HapMap pop.", ld.options$gts.source)
               else ld.options$gts.source,
            "// maxsnps =", ld.options$max.snps.per.window,
            "// rsquare >", ld.options$rsquare.min), 
     if(plot.variants) 
       paste("\nvariant options: vcf =",
             if(!is.null(var.options$vcf$file)) { 
               var.options$vcf$file
             } else {
               paste(var.options$vcf[[1]]$file, var.options$vcf[[2]]$file, sep = " ; ")
             },
             if(!is.null(var.options$vcf.info.prune)) paste("// filter", paste(var.options$vcf.info.prune, collapse = " ; ")), 
             if(!is.null(var.options$vcf.id.prune)) paste("// filter", paste(var.options$vcf.id.prune, collapse = " ; ")))
  )
  
  # region name (packet / panel name)
  # chr : startBP - endBP  (rsID)
  regions.df$packet.name <- factor(paste(
    regions.df$CHR, ":", 
    round(regions.df$region.startbp / 1000000, digits = 2), "M", "-", 
    round(regions.df$region.endbp / 1000000, digits = 2), "M", 
    "  (", regions.df$packet.snp, ")"
  ))
  
  # order of panels
  # we do not use a 'strip' function directly in the plot because we need the correct packet sorting
  # chromosome can be character or numeric -> make numbers equal length of the longest character length that occurs
  # e.g. chromosome 01 < 11
  chr.maxdigits <- max(nchar(as.vector(regions.df$CHR)))
  chr.print <- trim(as.character(regions.df$CHR))
  chr.print[nchar(chr.print) < chr.maxdigits] <- paste("0", chr.print[nchar(chr.print) < chr.maxdigits], sep = "")

  
  ##################### PLOT #####################


  plot.obj <- 
   xyplot(
    log10(P)*-1 ~ BP.mapped | packet.name,
    data = regions.df,
    group = pheno,
    # uuuhhh, this is the panel order by chr and bp - don't ask, just hope that it works...
    # unique keeps order of elements - first order regions.df by chr/pos, then get unique package names in the correct order
    index.cond = list(unique(regions.df$packet.name[order(as.vector(as.character(chr.print)), as.vector(as.numeric(regions.df$region.startbp)))])), 
    panel = function(...) {
      args <- list(...)
      packet.snp <- regions.df[args$subscripts[1], "packet.snp"]
      region.startbp <- regions.df[args$subscripts[1], "region.startbp"]
      region.endbp <- regions.df[args$subscripts[1], "region.endbp"]
      
      panel.xyplot(...)

      # manipulate with specific changes for this panel (e.g. smaller gene track size depending on number of elements)
      y.tracks.local <- y.tracks
      if(plot.genes) {
        gex <- genes.regions[[packet.snp]]
        if(!is.null(gex)) {
          if(gex$genes.density > 4) gex$genes.density <- 4
          gene.yextent <- (16 * ylim.space) * (1 + 3/gex$genes.density)
          panel.regionalplot.genelabels(
            gex$genes, 
            gex$genes.density, 
            exons = gex$exons, 
            region.startbp,
            region.endbp, 
            ystart = y.tracks.local["genes", "ystart"], 
            gene.yextent, 
            global.scale
          )
          ysize.adjusted <- gex$genes.density * gene.yextent
        } else {
          ysize.adjusted <- 0
        }
        # we have a dynamic y size of genes track
        # add excess space to the start position of the next track
        next.track.idx <- which(rownames(y.tracks.local) == "genes") +1
        y.tracks.local[next.track.idx, "ystart"] <- y.tracks.local["genes", "ystart"] - ysize.adjusted - y.tracks.local[next.track.idx, "yspace"]
      }
      
      if(plot.variants && !is.null(var.regions[[packet.snp]])) {
        panel.regionalplot.variants(
          var.regions[[packet.snp]], 
          x.start = region.startbp, 
          x.stop = region.endbp,
          y.start = y.tracks.local["var", "ystart"], 
          y.stop = y.tracks.local["var", "ystop"], 
          var.options, 
          out.format, 
          window.size,
          global.scale
        )
      }

      if(plot.ld) {
        panel.regionalplot.ld(
          snps = unique(regions.df[args$subscripts, c("SNP", "BP.mapped")]), 
          rsq = ld.regions[[packet.snp]], 
          ystart = y.tracks.local["ld", "ystart"], 
          ystop = y.tracks.local["ld", "ystop"], 
          global.scale, 
          rsquare.min = ld.options$rsquare.min, 
          show.rsquare.text = ld.options$show.rsquare.text
        )
      }
      
      if(is.data.frame(draw.snpname) || draw.snpname == "auto") {
        if(!is.data.frame(draw.snpname))
          # make auto-config
          draw.snpname <- data.frame(
            snps = snps[order(snps$BP, decreasing = T), "SNP"], 
            text = snps[order(snps$BP, decreasing = T), "SNP"],
            angle = if(nrow(snps) < 2) 60 else seq(20, 160, 140 / (nrow(snps) - 1)),
            length = 1, 
            cex = 1
          )
        snps.pos <- regions.df[args$subscripts, c("SNP", "BP.mapped", "P")]
        # for multiple pvalfiles, assign snpname to the graph with highest p-value
        snps.pos <- snps.pos[order(snps.pos$P), ]
        snps.pos <- snps.pos[!duplicated(snps.pos$SNP), ]
        snps.pos$P <- log10(snps.pos$P) * -1
        snps.pos <- merge(draw.snpname, snps.pos, by.x = "snps", by.y = "SNP")
        panel.regionalplot.snpnames(
          snps.pos, 
          out.format, 
          min(draw.snpname$cex) * global.scale
        )
      }

      panel.add(
        regions.df = regions.df, 
        packet.snp = packet.snp, 
        region.startbp = region.startbp, 
        region.endbp = region.endbp, 
        y.tracks = y.tracks, 
        global.scale = global.scale, 
        ylim.space = ylim.space,
        ...
      )
    }, 
    prepanel = function(x, y, subscripts, ...) { # define panel xsize (+ border)
      list(xlim = c(
        regions.df[subscripts[1], "region.startbp"] - window.size * 0.02, 
        regions.df[subscripts[1], "region.endbp"] + window.size * 0.02
      ))
    },
    ylim = c(ylim.lower - 20*ylim.space, ylim.upper + 30*ylim.space), 
    xlab = NULL, ylab = NULL, 
    between = list(x = 0, y = 1,5), 
    main = list(label = plot.title, cex = 0.8 * global.scale), 
    layout = c(1, out.format$panels.per.page), 
    as.table = TRUE, 
    lwd = 0.6 * global.scale,
    cex = 0.6 * global.scale, 
    type = "a", 
    par.strip.text = list(cex = 0.7 * global.scale, lineheight = 0.7 * global.scale),
    strip = function(bg, ...) strip.default(bg = "white", ...),
    key = {
      # custom key because we want small symbols and an extra legend when variants are used
      mykey <- simpleKey(levels(as.factor(regions.df$pheno)), cex = 0.6 * global.scale) 
      if(plot.variants && !is.null(var.options$vcf.info.colorize)) {
        mykey$text$lab <- c(mykey$text$lab, paste("variant:", var.options$vcf.info.colorize))
        mykey$points$col <- c(mykey$points$col, rainbow(length(var.options$vcf.info.colorize)))
        mykey$points$pch <- c(mykey$points$pch, rep(2,length(var.options$vcf.info.colorize)))
      }
      mykey$points$cex <- mykey$points$cex * 0.7 * global.scale
      mykey},
    xscale.components = function(...) {
      # scale baseposition on x axis ticks to MB
      dflt <- xscale.components.default(...)
      dflt$bottom$labels$labels <- as.numeric(dflt$bottom$labels$labels) / 1000000
      return(dflt)
    }, 
    scales = list(alternating = c(1), 
      x = list(relation = "free", tck = 0.25 * global.scale, cex = 0.5 * global.scale), 
      y = list(at = y.at, labels = y.labels, tck = 0.25 * global.scale, cex = 0.4 * global.scale)
    )
   )
  
  
  if(!is.null(out.format$file)) {
    message(paste("Writing plot to file", nextFilename("regionalplot", out.format$file), ""))
    if( out.format$file == "bmp" ) {
      trellis.device(
        bmp, 
        file = nextFilename("regionalplot", "bmp"), 
        width = out.format$paper.width, 
        height = out.format$paper.height, 
        units = "in", 
        res = 400, 
        pointsize = 1
        )
    } else {
      if( out.format$file == "jpeg" ) {
        trellis.device(
          jpeg, 
          file = nextFilename("regionalplot", "jpeg"), 
          width = out.format$paper.width, 
          height = out.format$paper.height, 
          units = "in", 
          res = 400, 
          pointsize = 1
          )
      } else {
        if( out.format$file == "png" ) {
          trellis.device(
            png, 
            file = nextFilename("regionalplot", "png"), 
            width = out.format$paper.width, 
            height = out.format$paper.height, 
            units = "in", 
            res = 400, 
            pointsize = 1
            )
        } else {
          trellis.device(
            pdf,
            file = nextFilename("regionalplot", "pdf"), 
            width = out.format$paper.width, 
            height = out.format$paper.height, 
            title = "Regional SNP association plot for multiple datasets"
            )
        }
      }
    }
    print(plot.obj)
    dev.off()
  }

  # cleanup : restore lattice options
  lattice.options(lattice.options.restore)

  return(plot.obj)

}

