% File man/postgwas.Rd
\name{postgwas}
\alias{postgwas}
\title{Default GWAS Visualization, Annotation and Network}

\description{
  Run a default set of postgwas functions for a quick analysis of your GWAS dataset.
}

\usage{
postgwas(
  gwas.dataset, 
  suggestive.p = 1e-06, 
  genomewide.p = 5e-08,
  gts.source = NULL, 
  biomart.config = biomartConfigs$hsapiens, 
  GOpackagename = "org.Hs.eg.db"
)
}

\arguments{
  \item{gwas.dataset}{A GWAS dataset. Several formats are accepted. Reads GenABEL, GEMMA, FAsTLMM or PLINK formats, or a default format containing columns 'SNP' and 'P' and optionally 'CHR' and 'BP'. The dataset can be passed as a filename (character argument) or a data frame containing the appropriate columns. For GenABEL, can also be a scan-gwaa object. Then, the P1df slot is used to provide p-values. In every case, each single dataset may not contain duplicate SNPs (is being checked). }
  \item{suggestive.p}{numeric(1). This p-value threshold selects SNPs that are considered for network analysis and regional plots.}
  \item{genomewide.p}{numeric(1). Currently, this only determines the height of the genomewide-line in the manhattanplot.}
  \item{gts.source}{vector(1). May contain genotypes of the GWAS study to allow calculation of LD information. Can be a HapMap population identifier (numeric) to retrieve genotyes for, or an object of class \code{\link{snp.data}} holding genotypes, or a GenABEL (.gwaa) or LINKAGE / PLINK (.ped) genotype file, with existing corresponding .phe and .map files, respectively. See the \code{gts.source} argument in the \code{\link{getGenotypes}} function for the exact file format specifications. }
  \item{biomart.config}{list. A configuration list for biomart data retrieval as described in \code{\link{biomartConfigs}}. Has to be changed when the species in the GWAS is not human.}
  \item{GOpackagename }{character(1). When GO enrichment in the network analysis should be used, this has to be the package name of an installed GO annotation package. It is preconfigured for human GO dataset, further species can be installed from \url{http://www.bioconductor.org/packages/2.10/data/annotation/} and the according package name passed to this argument when using non-human species. Can otherwise be NULL to avoid vertex colorization.}
}

\details{
 Applies the most important functions of postgwas sequentially to the supplied gwas dataset. This will produce several output files (see value section). Uses buffer variables (see \code{\link{postgwasBuffer}}), so all preexisting buffers are cleared. When the gts.source argument is not set, regionalplot lacks an LD track, and SNP to gene assignments are solely based on proximity. 
}

\value{
  Does not have a proper return value, but creates files as a side effect. 
  The following files are created / replaced: postgwas.snp2gene.csv, gwas2network.pdf, gwas2network.png, gwasGOenrich.csv, gwas2networkParams.txt, gwas2networkGraphVertices.csv, gwas2networkGraphEdges.csv
  The following files are consecutively numbered: regionalplot*.pdf, gwas2networkCommunities*.pdf, gwas2networkCommunities*.png
}

\seealso{
  \code{\link{snp2gene}}, \code{\link{regionalplot}}, \code{\link{gwas2network}}, \code{\link{removeNeighborSnps}}
}

\examples{
# this example requires online access 
# and the optional packages 'igraph' and 'GOSim' to be installed
\dontrun{
  gwas.dataset <- system.file("extdata", "example.gwas1", package = "postgwas")
  genotype.file <- system.file("extdata", "example.ped.xz", package = "postgwas")
  postgwas(
    gwas.dataset, 
    suggestive.p = 5e-02, 
    gts.source = genotype.file
  )
}
}