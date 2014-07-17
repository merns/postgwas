
# draws red LD triangle lines between SNPs in the y-range of the plot from y.start to y.stop
# Args:
#   snps:   Has to be a data frame with columns "SNP", "BP.mapped" containing snpIDs and base positions (x axis positions)
#   rsq:    Has to be a matrix containing r-square values. Columns and rows have to correspond to snps$SNP.
#   ystart: The y position of the LD triangles to start
#   ystop:  The y position of the LD triangles to end

panel.regionalplot.ld <- function(snps, rsq, ystart, ystop, global.scale, rsquare.min = 0.3, show.rsquare.text = FALSE) {

    if(is.null(rsq)) return()
  
    rsq[ rsq < rsquare.min ] <- NA

    # get pairs of snps where cells in rsq not NA
    ldsnps <- merge(data.frame(rownames(rsq)), data.frame(colnames(rsq)), all = TRUE)
    ldsnps <- ldsnps[!is.na(as.vector(rsq)), ]
    colnames(ldsnps) <- c("snp.x", "snp.y")

    if( nrow(ldsnps) > 0 ) {

      # set transparency (alpha) of to-draw lines according to r-square, normalized to [0,1]
      ldsnps$rsquare <- rsq[!is.na(as.vector(rsq))]
      ldsnps$alpha <- (ldsnps$rsquare - rsquare.min) / (max(ldsnps$rsquare) - rsquare.min + 0.01 )

      ldsnps <- merge(ldsnps, snps, by.x = "snp.x", by.y = "SNP")
      ldsnps <- merge(ldsnps, snps, by.x = "snp.y", by.y = "SNP")

      ldsnps$BP.mean <- ldsnps$BP.mapped.x - (ldsnps$BP.mapped.x - ldsnps$BP.mapped.y) /2
      ldsnps$BP.dist <- abs(ldsnps$BP.mapped.x - ldsnps$BP.mapped.y)

      # normalize triangle height to interval [ymin, ymax]
      ldsnps$height <- (ldsnps$BP.dist / max(ldsnps$BP.dist)) * (ystop-ystart) +ystart

      #draw triangles
      panel.segments( 
        x0 = t(ldsnps$BP.mapped.x), 
        y0 = ystart, 
        x1 = t(ldsnps$BP.mean), 
        y1 = t(ldsnps$height),
        col = "red", 
        alpha = ldsnps$alpha, 
        lwd = 0.7 * global.scale
      )
      panel.segments(
        x0 = t(ldsnps$BP.mapped.y), 
        y0 = ystart, 
        x1 = t(ldsnps$BP.mean), 
        y1 = t(ldsnps$height), 
        col = "red", 
        alpha = ldsnps$alpha, 
        lwd = 0.7 * global.scale
      )
      if( show.rsquare.text )
        panel.text( 
          x = t(ldsnps$BP.mean), 
          y = t(ldsnps$height), 
          labels = sub("^0", "", round(ldsnps$rsquare, digits = 2)),
          col = "red", 
          cex = 0.39 * global.scale
        )
    }
}
