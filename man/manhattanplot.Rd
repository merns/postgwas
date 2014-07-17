% File man/manhattanplot.Rd
\name{manhattanplot}
\alias{manhattanplot}
\title{Generate a Manhattan Plot from GWAS Data}

\description{
Produces a manhattan plot with highlighting and annotation features.
}

\usage{
manhattanplot(
  gwas.dataset, 
  highlight.logp = c(6, 7.3), 
  highlight.win = 125000, 
  highlight.color = NULL, 
  highlight.text = c("SNP", "genes"),
  highlight.cex = c(0.7, 1),
  highlight.fontface = c("italic", "bold"), 
  highlight.lines = FALSE, 
  ticks.y = FALSE,
  max.y = NULL,
  reduce.dataset = TRUE, 
  plot.title = NULL, 
  biomart.config = biomartConfigs$hsapiens, 
  use.buffer = FALSE, 
  toFile = nextFilename("manhattanplot", "pdf")
)
}

\arguments{
  \item{gwas.dataset}{A GWAS dataset. Several formats are accepted. Reads GenABEL, GEMMA, FAsTLMM or PLINK formats, or a default format containing columns 'SNP' and 'P' and optionally 'CHR' and 'BP'. The dataset can be passed as a filename (character argument) or a data frame containing the appropriate columns. For GenABEL, can also be a scan-gwaa object. Then, the P1df slot is used to provide p-values. For a PLINK file, only a single analysis model should be contained (e.g. one of DOM, REC, GENO).}
  \item{highlight.logp}{integer. The value defines a threshold for highlighting SNPs with -log(p) < highlight.logp. See also the arguments below. Can be a vector of multiple thresholds. }
  \item{highlight.win}{integer. For a SNP above the highlighting threshold, all SNPS lying within highlight.win basepair distance are also highlighted (applied to both sides). This argument is further used to determine a representative (i.e. lead) SNP of multiple proximate SNPs that all hit the highlighting threshold. The gene annotation is only plotted at the lead SNP of a block. The length of this argument should match the length of highlight.logp or is recycled if necessary.}
  \item{highlight.color}{character. Can be set to a user-defined color as listed in colors(). The length of this argument should match the length of highlight.logp or is recycled if necessary. }
  \item{highlight.text}{character. When 'SNP', draws SNP identifiers above the peak SNP of each highlighted locus. When "genes" runs snp2gene to annotate the closest genes in each direction (SNPs without gene annotation get their SNP ID plotted instead). When "none", draws no text annotation. The length of this argument should match the length of highlight.logp or is recycled if necessary.}
  \item{highlight.cex}{numeric. This is applied as cex argument (see \code{\link{par}}) to the gene labels, and in similar form to the plotted point symbol size of highlighted SNPs. The length of this argument should match the length of highlight.logp or is recycled if necessary.}
  \item{highlight.fontface}{character. Determines the fontface of the annotation text (can be "plain", "bold", "italic", "oblique", and "bold.italic"). The length of this argument should match the length of highlight.logp or is recycled if necessary. See details in the corresponding grid parameter documentation \code{\link{gpar}}. }
  \item{highlight.lines}{boolean. Determines whether a horizontal line is drawn at the \code{highlight.logp} thresholds. The length of this argument should match the length of highlight.logp or is recycled if necessary. }
  \item{ticks.y}{logical. When true, adds a graphical representation of the y axis to the plot (vertical line with tick marks).}
  \item{max.y}{numeric. The maximum -log(p) value to display. This can be useful when very small p's occur, cutting the data at a certain logp, or to extend the y axis to a larger range. When this argument is NULL, the y axis is scaled to the maximum / minimum p value from the gwas datafile.)}
  \item{reduce.dataset}{boolean(1). When TRUE, randomly removes SNPs to reduce data load, with increasing probablilty to be removed for low p-values. The top 0.2 percent of p-values remain untouched. Can also be a positive integer to do multiple removal cycles. }
  \item{plot.title}{character(1). Title text of the plot.}
  \item{biomart.config}{list. A biomart config list as described in \code{\link{biomartConfigs}}. This is only cosidered when highlight.text = "genes" and directly passed to the snp2gene function for gene annotations.}
  \item{use.buffer}{boolean(1). As in \code{\link{snp2gene}}.}
  \item{toFile}{character. A file/pathname to write the plot to (in pdf format). When NULL, plots to the current device (usually the display).}
}

\details{
Plots a grayscale-separated manhattan plot, with SNPs exceeding a user-defined threshold being highlighted, colorized and annotated with surrounding genes (or custom text). 
If a single locus gets genes annotated at several SNPs belonging to the same locus, extend the highlight.win size to unify them. 
}

\value{
 NULL, plots to the current device as a side effect.
}

\seealso{
  \code{\link{snp2gene}}
}

\examples{

\dontshow{
## offline usage: load and set buffer data for gene annotation
setPostgwasBuffer(
  snps = read.table(system.file("extdata", "postgwas.buffer.snps", package = "postgwas"), header = TRUE), 
  genes = read.table(system.file("extdata", "postgwas.buffer.genes", package = "postgwas"), header = TRUE)  
)
}
  
  # our file to plot
  gwas.dataset = system.file("extdata", "example.gwas1", package = "postgwas")
  
  # basic plot
  manhattanplot(
    gwas.dataset, 
    highlight.logp = c(3, 4), 
    use.buffer = TRUE, 
    plot.title = "truncated genomewide example dataset",
    toFile = NULL
  )
  dev.off()


\dontrun{
  # for yeast
  # we use a truncated synthetical dataset..
  gwas.dataset = system.file("extdata", "example.gwas3.xz", package = "postgwas")
  manhattanplot(
    gwas.dataset, 
    highlight.win = 500, 
    highlight.logp = c(3, 4), 
    biomart.config = biomartConfigs$scerevisiae, 
    use.buffer = FALSE
  )
  # the gene annotation overlaps a bit, either decrease the highlight.cex argument or 
  # change the plot width, as outlined in the next example  
}


  # we see that a locus appears twice (annotating the same genes)
  # increase window size to combine these loci
  # you can uncomment the 'tiff' and 'dev.off' to write a different format
  # the effects vary between operating systems:
  # experiment with the pointsize (= font size) and width/height (affects dots and symbols)
  # when the dots disappear the resolution is too low
\dontshow{
  if(capabilities("cairo")) {
    dev.new(type = "cairo", width = 12, height = 5)
  } else {
    dev.new(width = 12, height = 5)
  }
}
  # tiff("manhattanplot.tiff", width = 20, height = 8, units = "cm", res = 2000, pointsize = 10)
  manhattanplot(
    gwas.dataset, 
    highlight.logp = c(3, 4), 
    highlight.win = 130000,
    highlight.color = c("darkgreen", "red"),
    highlight.cex = c(0.7, 0.9),
    highlight.fontface = c("plain", "bold"),  
    use.buffer = TRUE,
    plot.title = "truncated genomewide example dataset", 
    toFile = NULL
  )
  # dev.off()
  
\dontshow{
  # cleanup
  clearPostgwasBuffer()
  rm(gwas.dataset)
}
}


