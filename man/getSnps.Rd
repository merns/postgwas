% File man/getSnps.Rd
\name{getSnps}
\alias{getSnpsByWin}
\alias{getSnpsByRS}
\title{Retrieve SNPs from local file or database with additional attributes}

\description{
Gets all SNPs within the specified region either from a mapfile or those available in the HapMap project. 
Gets all SNPs with position either from a mapfile or from the HapMap project. 
}

\usage{
getSnpsByWin(chr, bp.start, bp.end, map = NULL, pop.id = 2)
getSnpsByRS(snps, map = NULL, pop.id = 2)

}

\arguments{
  \item{chr}{character(1). A single chromosome number, formatted as in the referred mapfile or HapMap database}
  \item{bp.start}{numeric(1). Start base position of the window (inclusive)}
  \item{bp.end}{numeric(1). End base position of the window (inclusive)}
  \item{map}{When 
    \itemize{
    \item{NULL, will use hapmap snps }
    \item{character(1): fully qualified pathname, reads PLINK map file as in readMapfile.R }
    \item{data frame: assumes PLINK map file in map3 column format ("CHR", "SNP", "BP") }
  }
  }
  \item{pop.id}{numeric(1). HapMap population ID for which the SNPs are retrieved (should not have an effect but is mandatory, 2 = CEU)}
  \item{snps}{character. A vector of rs IDs.}
  \item{...}{May be parameters used by a subfunction (e.g. for LD calculation, can specify ped- and map parameters for genotype retrieval, see documentation of the function 'getGenotypes') }
}

\value{
byWin: A vector of SNPs
A vector of SNPs and columns "SNP", "CHR" and "BP"
}

\seealso{
  \code{\link{getGenotypes}}
}

\examples{
map.filename <- system.file("extdata", "example.map.xz", package = "postgwas")
getSnpsByWin(chr = 3, bp.start = 9300000, bp.end = 9500000, map = map.filename)
getSnpsByRS(snps = c("rs13240", "rs8411"), map = map.filename)
\dontshow{
  # cleanup
  rm(map.filename)
}
}
