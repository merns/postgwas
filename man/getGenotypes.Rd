% File man/getGenotypes.Rd
\name{getGenotypes}
\alias{getGenotypes}
\title{Retrieve Genotype Data}

\description{
Retrieves genotype data for given SNPs. Can do online HapMap queries or read ped/map (LINKAGE) or gwaa/phe (GenABEL) genotype files. For HapMap and LINKAGE mode, writes the retrieved genotype data to files as specified by the toFile argument.   
}

\usage{
getGenotypes(
  snps,
  gts.source = 2, 
  remove.homozygous = FALSE,
  return.characters = FALSE,
  return.genotypes = FALSE,
  toFile = "getGenotypes"
)
}

\arguments{
  \item{snps}{character. A vector with SNP IDs to get genotypes for.}
  \item{gts.source}{vector(1). This argument can be either a HapMap population ID (when numeric, e.g. 2 will download CEU sample genotypes), a genotype file in LINKAGE / PLINK format or a genotype file in GenABEL format (when type is character), or an object of class \code{\link{snp.data}}. For the GenABEL file format, either a .gwaa or .phe file can be specified, but both has to exist with the same basename in the same directory and will be read as specified by the \code{\link{load.gwaa.data}} function. For LINKAGE format, in the same way either a .ped or .map file can be specified, but both has to exist. The map file has to be in CHR, SNP, BP column format (PLINK option --map3). The pedigree (ped) file can have various formats, e.g. numeric or character nucleotides are accepted, space only delimiters or tab delimited with spaces between alleles etc. Additional columns ahead of genotype data (i.e. sex, family data, ...) can be included and will be ignored. LINKAGE files may also be compressed (see \code{\link{file}}), but in this case, both ped and map file have to be compressed. }
  \item{remove.homozygous}{boolean(1). Returns only genotypes of SNPs that have more than one allele state. }
  \item{return.characters}{boolean(1). When TRUE, the return value is always a matrix of character nucleotide values ACGT0 (or NULL for no data). }
  \item{return.genotypes}{boolean(1). When TRUE, the return value is always a list of \code{\link[genetics]{genotype}} objects (or NULL for no data). }
  \item{toFile}{character(1). Has no effect when gts.source is a GenABEL file, but has always be set for HapMap / LINKAGE file genotype retrieval. In the latter case, writes gwaa/phe and ped/map files of the retrieved genotype data, e.g. for later re-use. The value of this argument serves as basename of the files to be written. }
}

\value{
A \code{\link{snp.data}} object or NULL when genotypes are not available. May also return an error String when HapMap retrieval is selected and the connection is not available.  
Alternatively, when 
\itemize{
  \item{return.characters = TRUE: Returns a matrix of character values. Contains 2-character combinations of ACGT0 or NA values. SNPs are listed in columns.} 
  \item{return.genotypes = TRUE: Returns a list of \code{\link[genetics]{genotype}} objects. List elements are named by the SNP identifiers that have been retrieved.}
} 
}

\examples{
getGenotypes(
  snps = c("rs172154", "rs759704"), 
  gts.source = system.file("extdata", "example.gwaa", package = "postgwas"),
  return.characters = TRUE, 
  toFile = NULL
)
}
