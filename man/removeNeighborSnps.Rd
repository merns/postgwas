% File man/removeNeighborSnps.Rd
\name{removeNeighborSnps}
\alias{removeNeighborSnps}
\title{Remove SNPs within a Specified Range from a Data Frame}

\description{
This can be used to
\itemize{
  \item{select a representative SNP from a single chromosomal region}
  \item{ensure a certain distance between SNPs}
}
After application to a data frame of SNPs, this data frame will be pruned so that the remaining SNPs do not have a neighbor within the specified distance in the data frame. 
When the data frame contains a score column (p-value), the selected SNP(s) are asserted to have the lowest p-value (of the removed).
}

\usage{
removeNeighborSnps(
  snps, 
  maxdist = 50000, 
  biomart.config = biomartConfigs$hsapiens,
  use.buffer = FALSE
)
}

\arguments{
  \item{snps}{data frame. Has to contain the column 'SNP' containing a SNP identifier (mandatory) and optionally 'CHR', 'BP', 'P' for chromosome, base position and p-value (all numeric). When 'CHR' and 'BP' are not present, position information will be retrived from biomart using the specified config and / or buffer. The 'P' column is used as constraint for removal (see details) when present. May contain additional columns that are preserved.}
  \item{maxdist}{integer(1). Window size, SNPs with distance (<=) maxdist from another SNP are removed}
  \item{biomart.config}{list. Only needed when argument 'snps' lacks the columns 'CHR' and 'BP'. Contains values that are needed for biomaRt connection and data retrieval (i.e. dataset name, attribute and filter names for SNPs and genes). This is organism specific. Normally, a predefined list from \code{\link{biomartConfigs}} is specified, e.g. biomartConfigs$mmusculus for the mouse organism. Available predefined lists are shown with \code{names(biomartConfigs)}. See also \code{\link{biomartConfigs}}.}
  \item{use.buffer}{boolean(1). Only needed when argument 'snps' lacks the columns 'CHR' and 'BP'. When TRUE, buffers the downloaded snp annotation data in the \code{snps} buffer variable. When this variables already exists, data is not downloaded but used from there instead. See \code{\link{postgwasBuffer}} for more information on the buffer concept and variables. Make sure that pre-existing buffer data is current when use.buffer = TRUE. This facilitates the possibility to re-use data, alter data or provide custom data.}
}

\details{
This function removes from the data frame of snps all elements that have a neighbor within maxdist in this data frame. When the source data frame contains a column 'P', keeps always the SNP with the lowest p-value upon removal. Otherwise removes in order of rows in the data frame. Note that the results of the two methods may differ slightly (owing to moving the window by row or by p-value), though both are correct. For the p-value based version, the default maximum recursion depth of R might prohibit handling large datasets. Try setting options("expressions") = 50000 then. 
}

\value{
The source data frame 'snps', with neighbors removed.
}


\examples{
\dontrun{
  snps <- data.frame(SNP = c("rs2390012", "rs11161747"), P = c(0.04, 0.8))
  removeNeighborSnps(snps)
}
}


