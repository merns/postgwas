data.regionalplot.genes <- function(
  snps, 
  regions.df, 
  window.size, 
  gene.xspace, 
  geneds, 
  biomart.config, 
  use.buffer
) {

  message("Retrieving gene (and exon) information...")
  # get genes for all regions in a single rush
  genes.all <- bm.genes.regionalplot(
                    biomart.config, 
                    geneds, 
                    snps$CHR, 
                    snps$region.startbp, 
                    snps$region.endbp, 
                    use.buffer
                )
  genes.all <- na.omit(genes.all)
  genes.all <- genes.all[genes.all$name != "", ] # only interested in genes having a symbol

  # exons the same way
  if(exons.avail(biomart.config, use.buffer))
    exons.all <- bm.exons(
                  biomart.config, 
                  geneds, 
                  list(geneid = as.vector(genes.all$id), chr = genes.all$chromo), 
                  use.buffer
                )

  genes.regions <- lapply(
    1:nrow(snps),
    function(idx) {

      if(nrow(genes.all) < 1) return(NULL)
      genes <- genes.all[
                 genes.all$chromo == snps[idx, "CHR"] & 
                 # only genes that overlap with region (i.e. also partial matches)
                 (genes.all$start < as.numeric(snps[idx, "region.endbp"]) & 
                  genes.all$end > as.numeric(snps[idx, "region.startbp"]))
               , ]
      if(nrow(genes) < 1) return(NULL)

      if(exons.avail(biomart.config, use.buffer) && nrow(exons.all) > 0) {
        exons <- merge(exons.all, data.frame(chromo = genes$chromo, genestart = genes$start, geneend = genes$end))
      } else {
        exons <- NULL
      }

      # in the plot, each gene label should have half an inch of space (measured in bp)
      # additionally, when we have a global.scale factor, which influences character size, we have to apply it here as well
      # using gene.xspace, we calculate the maximum density of genes per xspace, so that the number of required y levels can be calculated in the panel
      density <- getDensity(genes, gene.xspace, "mean")
      ovlp <- merge(genes, genes, by = NULL) # cartesian product to check overlap between genes
      if(density < 2 & any(ovlp$start.x > ovlp$start.y & ovlp$start.x < ovlp$end.y)) density <- 2
      return(list(genes = genes, exons = exons, genes.density = density))

    }
  )
  names(genes.regions) <- snps$SNP
  return(genes.regions)
}
