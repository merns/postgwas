\name{postgwas-package}
\alias{postgwas-package}
\docType{package}
\title{
GWAS post-processing utilities
}

\description{
\itemize{
  \item{Annotate genes to SNPs (using proximity or LD information), see \code{\link{snp2gene}}}
  \item{Create \code{\link[=regionalplot]{regional}} and \code{\link[=manhattanplot]{manhattan}} plots}
  \item{Interaction network analysis of GWAS data, see \code{\link{gwas2network}}}
  \item{search literature (under development)}
}
}

\details{
This package contains several functions that are designed to simplify the interpretation of genomewide association studies (GWAS). 
Firstly, some frequently used visualization tools are implemented here in R that provide plots of publication-ready quality (manhattanplot, regionalplot). The manhattan plotting function offers highlighting features and the regional plots display genes, exons, p-values from multiple GWAS result files (comparative plot),linkage disequilibrium and SNP numbers upon users choice from diverse formats and data sources. There are options for additional rare variant track. Regional plots can be produced for many loci at once in a single pdf or png file. There is a function for comprehensive interaction network analysis of GWAS results, including download and construction of networks (gwas2network, getInteractions). 

All these functions can be run in a line with default parameters using the \code{\link{postgwas}} function, but it is recommended to use the individual functions for high-quality results.

Contributions to the package (additional functions, code extensions, bugfixed, suggestions for parameter changes etc) are much appreciated (please contact the package maintainer by mail). 
}

\author{
Milan Hiersche <mihi@uni-muenster.de>
}
