% File man/snp2gene.Rd
\name{snp2gene}
\alias{snp2gene}
\alias{snp2gene.prox}
\alias{snp2gene.LD}
\title{Annotate Genes to SNPs}

\description{
Adds columns with gene annotation to an input data frame of SNPs. Two modes are available: by proximity (genes covering the SNP and the next up- and downstream genes (including overlapping genes)) or, alternatively, genes that are in linkage equilibrium/disequilibrium with the query SNPs. Also adds biomart position annotation. 
}

\usage{
snp2gene.LD(
  snps, 
  gts.source = 2, 
  ld.win = 1000000, 
  biomart.config = biomartConfigs$hsapiens, 
  use.buffer = FALSE,
  by.genename = TRUE, 
  cores = 1
) 
snp2gene.prox(
  snps, 
  level = -1, 
  biomart.config = biomartConfigs$hsapiens, 
  use.buffer = FALSE,
  by.genename = TRUE,
  print.format = FALSE
)
snp2gene(
  snps, 
  gts.source = NULL, 
  biomart.config = biomartConfigs$hsapiens,
  use.buffer = FALSE, 
  by.genename = TRUE,
  level = -1, 
  ld.win = 1000000, 
  print.format = FALSE, 
  cores = 1
)
}

\arguments{
  \item{snps}{data frame. Has to contain a column named 'SNP' with rs-IDs (\url{http://www.ncbi.nlm.nih.gov/projects/SNP/}). May contain additional columns which are preserved. Will be made unique on the SNP column (duplicate SNPs discarded).}
  \item{gts.source}{vector(1). When NULL, annotates SNP by proximity. Otherwise, uses genotypes from a location as specified by the 'gts.source' argument to annotate genes by LD. Can take a numeric HapMap population ID to download genotypes for, or a GenABEL (.gwaa) or LINKAGE / PLINK format (.ped) genotype file (with corresponding .phe or .map file existing in the same directory), or an object of class \code{\link{snp.data}}. See the corresponding argument of the \code{\link{getGenotypes}} function for detailed specifications. GenABEL format is fastest and recommended. }
  \item{ld.win}{numeric(1). Specifies the total window size around a query SNP where genes are tested to be in LD with the query SNP. }
  \item{biomart.config}{list. Contains values that are needed for biomaRt connection and data retrieval (i.e. dataset name, attribute and filter names for SNPs and genes). This is organism specific. Normally, a predefined list from \code{\link{biomartConfigs}} is specified, e.g. biomartConfigs$mmusculus for the mouse organism. Available predefined lists are shown with \code{names(biomartConfigs)}. See also \code{\link{biomartConfigs}}.}
  \item{level}{numeric(1). When 0, returns for each SNP the (single) closest gene. When > 0, returns the level-closest gene in each direction, i.e. the closest in each direction for level = 1, and the second closest in each direction for level = 2 and so on. When < 0, returns the closest gene for each direction as for level = 1, plus all genes that overlap with the closest in each direction. }
  \item{print.format}{boolean(1). When FALSE (default), a decomposed data frame is returned, with one row per SNP-gene association. Contains an additional column for the direction of the SNP-gene relation (upstream, downstream, ...). LD annotation always uses the decomposed return format. When TRUE, returns the input data frame (one row per SNP) with three additional columns for covering, up- and downstream genes.}
  \item{cores}{integer(1). When > 1, parallelizes LD calculation and genotype retrieval by forking the parent R process. Number of concurrent processes is defined by the \code{cores} argument which has to be > 0. Might consume \code{cores}-fold memory. Parallelization only works on UNIX systems and might confuse the order of status messages.}
  \item{use.buffer}{boolean(1). When TRUE, buffers downloaded annotation data in the buffer variables \code{snps} and \code{genes}. When these variables already exists, data is not downloaded but used from there instead. This facilitates the possibility to re-use data, alter data or provide custom data. See \code{\link{postgwasBuffer}} for more information on the buffer concept and variables. When setting use.buffer = TRUE, make sure that existing buffer data is current!}
  \item{by.genename}{boolean(1). When TRUE, only genes that have a name are annotated. For genes with a single name but multiple IDs, a random representative ID is selected. When FALSE, uses all genes.}
}

\details{
  When gene annotation by proximity is requested, we have to deal with cases of overlapping genes. This can be many (e.g. TRBV gene variants or dense gene clusters). 
  When level < 0, all overlapping genes will be annotated (capped by cluster.maxsize) and form multiple rows in the decomposed format or, for the print format, occur as single string where gene names are separated by '|'. 
  LD calculation uses several measures to assess the amount of linkage disequilibrium between a SNP and a gene. In a preliminary step, the r-square correlation between the queried SNP and SNPs within the gene is calculated (uses the \code{\link{r2fast}} function of GenABEL). 
  Due to efficency constraints, only genes within a user-defined window (default: 1 MB) are considered (genes that overlap partially with the window are included, and the nearest gene is always considered regardless of the window parameter), and there is a cap of 100 SNPs per gene for LD calculation (when there are more, 100 SNPs are selected by position, evenly distributed over the gene). Output then contains 
  \itemize{
    \item{the best r-square value between the query SNP and all SNPs from the gene }
    \item{the mean r-square value between the query SNP and all SNPs from the gene }
    \item{the standard deviation of r-square value between the query SNP and all SNPs from the gene }
  }
  When the \code{gts.source} is a HapMap population or linkage file, all genotypes retrieved will be written to the files snp2gene.gwaa, snp2gene.phe, snp2gene.ped and snp2gene.map in the current directory.  
The whole annotation process makes use of the biomaRt database interface. All positions used and returned refer to the genome build version as provided by the according biomaRt retrival. The database and variables used can be specified by the user, e.g. it is possible to annotate accession numbers by altering the biomart.config, setting \code{myconfig$gene$attr$name <- 'uniprot_swissprot_accession'} 
}

\value{
  A data frame with all columns from the query data frame 'snps' plus columns 'BP' for the biomart position (renames an existing BP column to BP.original) and
  \itemize{
    \item{for annotation by proximity, columns 'geneid', 'genename', 'start', 'end' containing the gene id/name and position as in the biomart config and 'direction' giving the mode of annotation (cover, up- or downstream) }
    \item{for annotation by LD, columns 'geneid', 'genename', 'start', 'end' containing the gene id/name and position as in the biomart config, and 'ld.mean' and 'ld.max' with the according rsquare value}
  }
  Each association of a SNP to a gene forms a row of the data frame. 
}

\seealso{
  \code{\link{biomartConfigs}}
}

\examples{

snps <- data.frame(SNP = c("rs172154", "rs759704"))

\dontshow{
## offline usage: load and set buffer data
setPostgwasBuffer(
  snps = read.table(system.file("extdata", "postgwas.buffer.snps", package = "postgwas"), header = TRUE), 
  genes = read.table(system.file("extdata", "postgwas.buffer.genes", package = "postgwas"), header = TRUE)  
)
## offline LD annotation needs genotype files, can be ped/map or gwaa/phe
genotype.file <- system.file("extdata", "example.gwaa", package = "postgwas")
}

\dontrun{
## simplest annotation, online
annot.prox <- snp2gene.prox(snps)
annot.ld <- snp2gene.LD(snps)

## extract names of genes in LD
geneset <- annot.ld[annot.ld$ld.max > 0.3, "genename"]

## the same with a different organism
snps <- data.frame(SNP = c("s01-10027", "s04-1469331", "s10-240843", "s15-474479"))
snp2gene.prox(snps, biomart.config = biomartConfigs$scerevisiae)

## how to use protein accession numbers
bconf <- biomartConfigs$hsapiens
bconf$gene$attr$name <- "uniprot_swissprot_accession"
snp2gene.prox(snps, biomart.config = bconf)

## this shows a list of further biomart attributes for human genes:
listAttributes(useMart("ensembl", "hsapiens_gene_ensembl"))
}

## annotation using buffer variables (runs offline when set properly)
snp2gene.LD(snps, gts.source = genotype.file, use.buffer = TRUE)
snp2gene.prox(snps, use.buffer = TRUE)

\dontshow{
  # cleanup
  clearPostgwasBuffer()
  rm(genotype.file, snps)
}

}

