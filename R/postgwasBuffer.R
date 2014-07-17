# the environment postgwasBuffer is constructed in the onLoad function. 
# all buffer variables in the environment postgwasBuffer always have to exist (and do so by implementation)
# they can be NULL which means 'not set'

clearPostgwasBuffer <- function() {
  assign("snps", NULL, envir = postgwasBuffer)
  assign("genes", NULL, envir = postgwasBuffer)
  assign("genes.regionalplot", NULL, envir = postgwasBuffer)
  assign("exons.regionalplot", NULL, envir = postgwasBuffer)
  assign("ld.regionalplot", NULL, envir = postgwasBuffer)
  assign("goterms", NULL, envir = postgwasBuffer)
}


getPostgwasBuffer <- function() {
  return(list(
          snps = get("snps", envir = postgwasBuffer), 
          genes = get("genes", envir = postgwasBuffer), 
          genes.regionalplot = get("genes.regionalplot", envir = postgwasBuffer), 
          exons.regionalplot = get("exons.regionalplot", envir = postgwasBuffer), 
          ld.regionalplot = get("ld.regionalplot", envir = postgwasBuffer), 
          goterms = get("goterms", envir = postgwasBuffer)
      ))
}


setPostgwasBuffer <- function(
    uselist = FALSE, 
    snps, 
    genes, 
    genes.regionalplot, 
    exons.regionalplot, 
    ld.regionalplot, 
    goterms
) {
  
  if(is.list(uselist)) {
    if(!all(names(uselist) %in% c("snps", "genes", "genes.regionalplot", "exons.regionalplot", "ld.regionalplot", "goterms"))) 
      stop("setPostgwasBuffer: The supplied list of buffer variables may only contain elements named snps, genes, genes.regionalplot, exons.regionalplot, ld.regionalplot, goterms")
    if(!is.null(uselist$snps))
      snps <- uselist$snps
    if(!is.null(uselist$genes))
      genes <- uselist$genes
    if(!is.null(uselist$genes.regionalplot))
      genes.regionalplot <- uselist$genes.regionalplot
    if(!is.null(uselist$exons.regionalplot))
      exons.regionalplot <- uselist$exons.regionalplot
    if(!is.null(uselist$ld.regionalplot))
      ld.regionalplot <- uselist$ld.regionalplot
    if(!is.null(uselist$goterms))
      goterms <- uselist$goterms
  } 
  
  if(!missing(snps))
    if(!is.data.frame(snps) && !is.null(snps))
      stop("setPostgwasBuffer: Wrong format of argument 'snps'")
  if(!missing(genes))
    if(!is.data.frame(genes) && !is.null(genes))
      stop("setPostgwasBuffer: Wrong format of argument 'genes'")
  if(!missing(genes.regionalplot))
    if(!is.data.frame(genes.regionalplot) && !is.null(genes.regionalplot))
      stop("setPostgwasBuffer: Wrong format of argument 'genes.regionalplot'")
  if(!missing(exons.regionalplot))
    if(!is.data.frame(exons.regionalplot) && !is.null(exons.regionalplot))
      stop("setPostgwasBuffer: Wrong format of argument 'exons.regionalplot'")
  if(!missing(ld.regionalplot))
    if(!is.list(ld.regionalplot) && !is.null(ld.regionalplot))
      stop("setPostgwasBuffer: Wrong format of argument 'ld.regionalplot'")
  if(!missing(goterms))
    if(!is.data.frame(goterms) && !is.null(goterms))
      stop("setPostgwasBuffer: Wrong format of argument 'goterms'")
  
  if(!missing(snps))
    assign("snps", snps, envir = postgwasBuffer)
  if(!missing(genes))
    assign("genes", genes, envir = postgwasBuffer)
  if(!missing(genes.regionalplot))
    assign("genes.regionalplot", genes.regionalplot, envir = postgwasBuffer)
  if(!missing(exons.regionalplot))
    assign("exons.regionalplot", exons.regionalplot, envir = postgwasBuffer)
  if(!missing(ld.regionalplot))
    assign("ld.regionalplot", ld.regionalplot, envir = postgwasBuffer)
  if(!missing(goterms))
    assign("goterms", goterms, envir = postgwasBuffer)
}
