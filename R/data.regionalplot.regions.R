data.regionalplot.regions <- function(snps, p.all, window.size) {

  # center regions around query snp
  # we have assured that query snp has position from pval file
  message("Assembling regions...")
  regions <- sapply(
    as.vector(snps$SNP), 
    function(snp) {
      find.chr <- as.vector(snps[snps$SNP == snp, "CHR"])
      find.bp <- as.numeric(as.vector(snps[snps$SNP == snp, "BP"]))
      region <- p.all[
        p.all$CHR == find.chr & 
        as.numeric(as.vector(p.all$BP)) < find.bp + window.size/2 & 
        as.numeric(as.vector(p.all$BP)) > find.bp - window.size/2
      , ]
      region$SNP <- as.vector(region$SNP) # remove spurious levels from data frame
      region$packet.snp <- snp
      return(region)
    }, 
    simplify = F
  )
  regions.df <- list2df(regions)
  
}