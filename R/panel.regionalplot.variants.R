variants.hist <- function(
                    variants, 
                    y.start, 
                    y.stop, 
                    af.max, 
                    dual.vcf, 
                    scale.factor, 
                    alpha = 1
                  ) {

  if(nrow(variants) < 1) return()

  if(dual.vcf) {
    # comparative histogram between two datasets
    # but only for a diff > 0, otherwise draw horizontal line
    variants.diff <- variants[abs(variants$AF.x - variants$AF.y) > 0.01, ]
    if(nrow(variants.diff) > 0)
      panel.arrows(
        x0 = variants.diff$BP.mapped, 
        x1 = variants.diff$BP.mapped,
        y0 = y.stop + (y.start - y.stop) * (variants.diff$AF.x / af.max),
        y1 = y.stop + (y.start - y.stop) * (variants.diff$AF.y / af.max), 
        col = variants.diff$COLOR, 
        alpha = alpha, 
        length = 0.02 * scale.factor, 
        lwd = 1 * sqrt(scale.factor)
      )
    variants <- variants[abs(variants$AF.x - variants$AF.y) <= 0.01, ]
    panel.points(
      x = variants$BP.mapped,
      y = y.stop + (y.start - y.stop) * (variants$AF.x / af.max),
      col = variants$COLOR,
      pch = "-",
      alpha = alpha,
      cex = 0.5 * scale.factor
    )
  } else {
    # draw standard histogram
    panel.segments(
      x0 = variants$BP.mapped, 
      x1 = variants$BP.mapped,
      y0 = y.stop,
      y1 = y.stop + (y.start - y.stop) * (variants$AF / af.max), 
      col = variants$COLOR, 
      alpha = alpha, 
      lwd = 1 * sqrt(scale.factor)
    )
    panel.points(
      x = variants$BP.mapped,
      y = y.stop + (y.start - y.stop) * (variants$AF / af.max),
      col = variants$COLOR,
      pch = "-", 
      alpha = alpha,
      cex = 1 * scale.factor
    )
  }
}

variants.line <- function(variants, y.stop, scale.factor) {
  # add thin vertical calibration lines across panel
  if(nrow(variants) == 0) return()
  panel.segments(
    x0 = variants$BP.mapped, 
    x1 = variants$BP.mapped,
    y0 = 100,    # start somewhere at the very top of the plot
    y1 = y.stop, 
    col = variants$COLOR,
    lwd = 0.6 * sqrt(scale.factor), 
    lty = 3
  )
}

variants.labels <- function(variants, dual.vcf, y.start, y.stop, window.size, out.format, scale.factor) {
  if(nrow(variants) == 0) return()
  label.xspace <- (window.size / out.format$paper.width) * 2 * scale.factor  # two inches per label
  dens <- getDensity(variants, label.xspace, "BP.mapped")
  # we scale the label size quadratically, so we can lower the needed space accordingly
  if(dens > 4)
    dens <- getDensity(variants, label.xspace / (dens)^(1/3), "BP.mapped")
  total.space <- y.stop - y.start
  space.per.label <- total.space / dens
  offset <- space.per.label / 2

  variants$ylevel <- rep(0:(dens-1) * space.per.label + y.start, nrow(variants))[1:nrow(variants)]

  panel.text(
    x = variants$BP.mapped, 
    y = variants$ylevel + offset, 
    labels = paste(
               variants$ID, 
               variants$BP.mapped, 
               if(dual.vcf) paste(variants$BP.x, variants$BP.y, sep = " : ") else variants$BP, 
               sep = " : "
              ),
    cex = (if(dens <= 3) 0.35  else 1.15/dens) * scale.factor^2, 
    col = variants$COLOR, 
    pos = 4
  )
}

panel.regionalplot.variants <- function( 
                    variants, # data.frame(ID, BP.mapped, BP, COLOR, AF) ; AF can also be AF.x and AF.y
                    x.start, 
                    x.stop, 
                    y.start,
                    y.stop,
                    var.options, 
                    out.format, 
                    window.size, 
                    scale.factor
                  ) {

  dual.vcf <- if("AF.y" %in% colnames(variants)) TRUE else FALSE
  variants.colored <- variants[variants$COLOR != "#000000", ]
  variants.black <- variants[variants$COLOR == "#000000", ]
  y.divisor <- y.stop + (y.start - y.stop) / 3
  
  # the lower part of the plot has a space separated by horizontal lines where the vertical variant lines end
  # the end position within this space determines the allele count (histogram-like)
  # even below, the variant identifier and allele count are drawn

  panel.segments(
    x0 = x.start, 
    x1 = x.stop, 
    y0 = c(y.start, y.stop), 
    y1 = c(y.start, y.stop), 
    col = "black", lty = 2, lwd = 0.5 * scale.factor
  )
  
  if(var.options$details <= 2) {
    variants.hist(variants.black, y.start, y.stop, var.options$vcf.af.prune, dual.vcf, scale.factor, alpha = var.options$hist.alpha)
    variants.hist(variants.colored, y.start, y.stop, var.options$vcf.af.prune, dual.vcf, scale.factor)
  }

  if(var.options$details >= 2) {
    variants.line(
      variants.colored, 
      y.stop,
      scale.factor
    )
  }
  
  if(var.options$details >= 3) {
    # break available y space into histogram and label parts
    panel.segments(
      x0 = x.start, 
      x1 = x.stop, 
      y0 = y.divisor, 
      y1 = y.divisor, 
      col = "black", alpha = 0.2, lty = 2, lwd = 0.5 * scale.factor
    )
    variants.hist(
      variants.black,
      y.start = y.divisor, 
      y.stop = y.stop, 
      var.options$vcf.af.prune,
      dual.vcf,
      scale.factor,
      alpha = var.options$hist.alpha
    )
    variants.hist(
      variants.colored, 
      y.start = y.divisor, 
      y.stop = y.stop, 
      var.options$vcf.af.prune,
      dual.vcf,
      scale.factor
    )
    variants.labels(
      variants.colored,
      dual.vcf, 
      y.start = y.start, 
      y.stop = y.divisor,
      window.size,
      out.format,
      scale.factor
    )
  }
  
}

