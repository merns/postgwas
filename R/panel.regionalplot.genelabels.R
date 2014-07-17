
# Shows an arrow (orientation = strand) and the according gene symbol in the plot.
# Can also mark exons with a thickening of the gene line.
#
# Args:
#   genes:        a data frame of genes to plot, with columns "start", "end", "mean", "name", "strand" ; strand = -1 means reverse, everything else forward 
#   minbp:        an integer number defining the smallest base position of the plot packet (~xlim)
#   minbp:        an integer number defining the largest base position of the plot packet (~xlim)
#   plim:         y extent
#   line.ypos:    y position of the line marking the gene
#   labels.ypos:  y position of the gene label (first character)
#   exons:        a data frame containing columns "exon_chrom_start" and "exon_chrom_end" with start and end positions of exons
#   scale.factor: character expansion value as in e.g 'xyplot'. Is also used for other scaling calculations (e.g. exon size, lwd) here.
#   truncate:     When TRUE, gene start / end positions will be truncated to xstart and xstop. Exons with exon_chrom_start > xstop
#                  or exon_chrom_end < xstart are dicarded (that are those that lie completely beyond xstart and xstop).

panel.regionalplot.genelabels <- function(
                      genes, 
                      levels, 
                      exons = NULL, 
                      xstart, 
                      xstop, 
                      ystart, 
                      ygenesize, 
                      scale.factor, 
                      truncate = TRUE,
                      ...
                    ) {

  if(nrow(genes) > 0) {

    genes <- genes[order(genes$mean), ]
    genes$ypos <- rep(-1 * (ygenesize/3.5 + ygenesize * 0:(levels-1)) + ystart, nrow(genes))[1:nrow(genes)]

    # transfer y position to exon
    if(!is.null(exons) && nrow(exons) > 0)
      exons <- merge(exons, genes, by.x = c("chromo", "genestart", "geneend"), by.y = c("chromo", "start", "end"))
    
    if(truncate) {
      if(!is.null(exons) && nrow(exons) > 0)
        exons <- exons[exons$end > xstart & exons$start < xstop, ]
      genes[genes$start < xstart, "start"] <- xstart
      genes[genes$end > xstop, "end"] <- xstop
      # do not show gene name when mean is outside window
      genes[genes$mean < xstart | genes$mean > xstop, "name"] <- ""
    }
	
    # gene arrows - forward and reverse
    genes.rev <- genes[genes$strand == -1, ]
    if(nrow(genes.rev) > 0)
      panel.arrows( 
        x0 = genes.rev$end, 
        y0 = genes.rev$ypos, 
        x1 = genes.rev$start, 
        y1 = genes.rev$ypos,
        col = "forestgreen",
        length = 0.05 * scale.factor^2,
        lwd = scale.factor, 
        ...
      )

    genes.fwd <- genes[genes$strand != -1, ]
    if(nrow(genes.fwd) > 0)
      panel.arrows( 
        x0 = genes.fwd$start, 
        y0 = genes.fwd$ypos, 
        x1 = genes.fwd$end, 
        y1 = genes.fwd$ypos, 
        col = "forestgreen",
        length = 0.05 * scale.factor^2,
        lwd = scale.factor,
        ...
      )

    # y position of exons == gene arrow positions
    if(!is.null(exons) && nrow(exons) > 0) {
      panel.rect( 
        xleft = exons$start,
        ybottom = exons$ypos - 0.2 / scale.factor^2,
        xright = exons$end,
        ytop = exons$ypos + 0.2 / scale.factor^2,
        border = NA, 
        col = "black"
      )
    }


    # gene names
    panel.text( 
      x = genes$mean, 
      y = genes$ypos, 
      labels = genes$name,
      col = "forestgreen",
      cex = 0.4 * scale.factor, 
      pos = 1,
      ...
    )

  }
	
}

