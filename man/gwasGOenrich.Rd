% File man/gwasGOenrich.Rd
\name{gwasGOenrich}
\alias{gwasGOenrich}
\title{Calculate Enrichment of GO Terms in GWAS Datasets}

\description{
This function wraps functionality from the topGO package, providing a simplified interface for GO term enrichment for GWAS datasets. Before running this function on a GWAS dataset, a SNP to gene mapping and calculation of gene-wise p-values has to be applied to the GWAS dataset, e.g. as offered by the \code{\link{snp2gene}} and \code{\link{gene2p}} functions. 
}

\usage{
gwasGOenrich(
    gwas, 
    ontology = "BP", 
    pkgname.GO = "org.Hs.eg.db", 
    topGOalgorithm = "weight01",
    pruneTermsBySize = 10, 
    plotSigTermsToFile = 5
)
}

\arguments{
  \item{gwas}{data frame. Has to contain a column of gene-wise p-values (see also \code{\link{gene2p}} named 'gene.p' and a column named 'geneid' or 'genename' containing entrez IDs, ensembl gene IDs, gene symbols or names, respectively.}
  \item{ontology}{character(1). Can be "BP", "MF" or "CC" to use the according ontology for enrichment (i.e. biological process, molecular function, cellular component).}
  \item{pkgname.GO}{character(1). Name of the GO annotation package to use. This is 'org.Hs.eg.db' for human datasets, but alternative 'org...' packages which are available via bioconductor for different organisms can be specified.}
  \item{topGOalgorithm}{character(1). Is passed to the \code{algorithm} argument of topGO. Possible values are listed by running \code{\link[topGO]{whichAlgorithms}}. E.g. 'elim' decorrelates the GO tree so that the most specific out of several similar terms is contained in the results, and 'weight01' has a similar effect. The 'classic' option will just run the test on all terms. }
  \item{pruneTermsBySize}{numeric(1). Is passed to the \code{nodeSize} argument of topGO. When >1, all GO terms with less than pruneTermsBySize genes will be discarded before the enrichment.}
  \item{plotSigTermsToFile}{numeric(1). When > 0, the topGO method \code{printGraph} will be called to produce a file containing a graph plot of the most significant terms and its descendant terms, together with enrichment stats. The storeGraphPlot most significant terms will be contained in the plot.}
}

\details{
This function builds a topGO data object at hand of the supplied gwas dataset and applies \code{\link[topGO]{runTest}} with the \code{statistic = "ks"} option which runs a gene-set enrichment like analyis. 
}

\value{
A data frame with columns 'go_id', 'Term' and 'P' containing identifier, description and enrichment p-value of all GO terms tested.
When \code{plotSigTerms} is > 0, a pdf file (name starting with 'enrichment...) will be produced as a side effect containing the graph plot for selected terms produced by topGO.  
}


\references{
  Adrian Alexa and Jorg Rahnenfuhrer (2010). topGO: topGO: Enrichment analysis for Gene Ontology. R package version 2.8.0.
}

