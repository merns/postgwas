% File man/bm.remapSNPs.Rd
\name{bm.remapSnps}
\alias{bm.remapSnps}
\title{Annotate Position to SNPs}

\description{
  Annotate current biomart base position information to a file containing SNPs.
}

\usage{
bm.remapSnps(
  file, 
  snp.col = "SNP", 
  config = biomartConfigs$hsapiens, 
  ds = bm.init.snps(config), 
  use.buffer = FALSE, 
  toFile = paste(file, ".remapped", sep = "")
)
}

\arguments{
  \item{file}{character(1). Path/filename, has to be readable by read.table (default settings). Has to contain a column of SNP identifiers that is used for annotation (values should match the selected biomart SNP identifier field, i.e. rs-IDs). See the 'snp.col' argument on how this column is identified. Rows that contain NA values are removed.}
  \item{snp.col}{character(1) or numeric(1). Identifies the column containing SNP identifiers in the source file. When this is a numeric argument, the source file will be read assuming no header line, using the snp.col-th column for SNP identifiers (fist column = 1). When snp.col is a character value, a header line is assumed in the source file, and the column name matching snp.col is taken as containing SNP identifiers.}
  \item{config}{A list containing values that are needed for biomaRt connection and data retrieval. See \code{\link{biomartConfigs}}.}
  \item{ds}{Optionally, a biomart connection object for re-use scenarios to avoid overhead creating connections.}
  \item{use.buffer}{Uses the buffer variable \code{snps} for annotation when it exists and contain valid data in biomart format. When it does not exists, writes the biomart data used to that variable. See \code{\link{postgwasBuffer}} for more information on buffer variables. }
  \item{toFile}{character(1). Path/filename, destination to write the remapped data to. Disabled when NULL.}
}

\value{
  A data frame containing the original file data plus columns 'CHR' and 'BP' with the biomart position annotation. Existing 'BP and CHR' columns will be renamed with the suffix .original. Unknown SNPs will be preserved and contain NA in the added columns. Order of SNPs is being preserved.  
}

\seealso{
  \code{\link{biomartConfigs}}
}
