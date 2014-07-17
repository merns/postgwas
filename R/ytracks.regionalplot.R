ytracks.regionalplot <- function(
                          ylim.upper, 
                          ylim.space, 
                          plot.genes = TRUE, 
                          plot.ld = TRUE, 
                          plot.variants = FALSE, 
                          var.options = NULL
                        ) {
  # return a data frame of track info with columns ysize, yspace, ystart, ystop
  # the ysize is the available space for the features to be drawn, 
  # and yspace accounts for empty space preceding the data track und thus allowing a visual separation on tracks
  # both are positive values (though accuring on the negative scale)
  # ystart and ystop are the (negative valued) positions on the y axis for the track, (without yspace)
  # rownames should correspond to the upcoming y axis labels

  tracks <- data.frame(
    ysize =  0, 
    yspace = 10 * ylim.space, 
    ystart = -10 * ylim.space, 
    ystop =  -10 * ylim.space, 
    name = "feature tracks:", 
    row.names = "default", 
    stringsAsFactors = F
  )
  
  if(plot.genes) {
    tracks <- rbind(tracks, "genes" = list(
      ysize  =  ylim.upper, 
      yspace =  10 * ylim.space, 
      ystart = -10 * ylim.space, 
      ystop  = -(ylim.upper + 10 * ylim.space), 
      name = "genes"
    ))
  }

  if(plot.ld) {
    tracks <- rbind(tracks, "ld" = list(
      ysize  =  0.6 * ylim.upper, 
      yspace =  15 * ylim.space, 
      ystart =  -(sum(tracks$ysize, tracks$yspace) + 15 * ylim.space), 
      ystop  =  -(sum(tracks$ysize, tracks$yspace) + 15 * ylim.space + 0.6 * ylim.upper), 
      name = "LD (r\u00B2)"
    ))
  }
  
  if(plot.variants) {
    ysize <- if(var.options$details <= 2) 0.5 * ylim.upper else ylim.upper
    tracks <- rbind(tracks, "var" = list(
      ysize  =  ysize, 
      yspace =  15 * ylim.space, 
      ystart =  -(sum(tracks$ysize, tracks$yspace) + 15 * ylim.space),
      ystop  =  -(sum(tracks$ysize, tracks$yspace) + 15 * ylim.space + ysize),
      name = paste("variants (0-", var.options$vcf.af.prune, "% AF)", sep = "")
    ))
  }

  return(tracks)
}
