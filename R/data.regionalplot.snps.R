data.regionalplot.snps <- function(snps, p.all) {

# require SNPs to occur in pval file(s)
  if( !all(snps$SNP %in% p.all$SNP) )
    stop(paste(
      "The following SNPs do not occur in the supplied p-value lists:", 
      paste(snps[!snps$SNP %in% p.all$SNP, "SNP"], collapse = ",")
   ))

  # add position from pval files to snps  
  snps <- p.all[p.all$SNP %in% as.vector(snps$SNP), 1:3]
  # remove spurious levels from data frame
  snps$SNP <- as.vector(snps$SNP)

  # make SNPs unique 
  snps <- snps[!duplicated(snps$SNP), ]
  
  return(snps)
  
}
