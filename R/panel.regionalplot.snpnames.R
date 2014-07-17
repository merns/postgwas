
# Draw a small connector line with the snp name at the correspongin position in the graph.
# An angle and length parameter can be specified for every SNP that controls the connector orientation and length.
#
# Args:
#   snps:   a data frame with columns "snps", "BP", "P", "angle", "length", "cex" 
#   x.size: the x axis size (units as used in the panel) 
#   y.size: the total y axis heigth (units as used in the panel) 
#   ratio:  the ratio of x to y panel dimensions (in actual pixels)
# 	cex:    scaling parameter for label size and connector thickness

panel.regionalplot.snpnames <- function(snps, out.format, cex) {

  x.size <- abs(diff(current.panel.limits()$xlim))
  y.size <- abs(diff(current.panel.limits()$ylim))
  ratio <- (out.format$paper.width / 0.82) / (out.format$paper.height / out.format$panels.per.page)
  
	snps$angle <- abs(snps$angle %% 360)
	# connector length = a 15th fraction of panel y size
	connector.ylength <- (y.size / 15) * snps$length
	connector.xlength <- (x.size / (15 * ratio)) * snps$length
	snps$angle.rad <- pi/180 * snps$angle
	
	# apply the rotation to the connector: sin is the y axis offset, cos the x axis
	# also truncate the connector length a bit to leave some space between connector and graph (and connector an label)
	panel.segments(
			snps$BP + (connector.xlength / 5) * cos(snps$angle.rad), 
			snps$P - (connector.ylength / 5) * -sin(snps$angle.rad), 
			snps$BP + (connector.xlength - connector.xlength / 8) * cos(snps$angle.rad), 
			snps$P - (connector.ylength - connector.ylength / 8) * -sin(snps$angle.rad), 
			lwd = 0.7 * cex
	)
	
	# discriminate snp labels by coordinate system quadrants they lie in
	# use different label offsets for each quadrant
	lapply(
			1:4, 
			function(quadr) {
				snps.quadrant <- snps[ snps$angle >= 90 * (quadr-1) & snps$angle < 90 * quadr, ]
				if(quadr == 1)
					adj <- c(0, 0)
				if(quadr == 2)
					adj <- c(1, 0)
				if(quadr == 3)
					adj <- c(1, 1)					
				if(quadr == 4)
					adj <- c(0, 1)
				panel.text(
						snps.quadrant$BP + connector.xlength * cos(snps.quadrant$angle.rad), 
						snps.quadrant$P - connector.ylength * -sin(snps.quadrant$angle.rad), 
						as.vector(snps.quadrant$text), 
						adj = adj, 
						font = 3, # italic 
						cex = 0.5 * cex
				)
			}
	)

}
