% File man/gene2p.Rd
\name{gene2p}
\alias{gene2p}
\alias{GATES}
\alias{SpD}
\title{Aggregate P-Values for Genes with multiple SNPs}

\description{
Calculates aggregate (combined) p-values for genes assigned to multiple SNPs, by taking into account 
the dependency structure between SNPs (mainly based on LD). This is done using the wrapping function 
gene2p, which can take specific algorithm implementations as function argument. Currently, re-implementations 
of the GATES and SpD algorithms (see references for original publications and authors) are available. The 
wrapping function 'gene2p' will be called by the user and manages LD calculation (\code{\link{r2fast}} function 
of GenABEL) and parallelization with subsequent application of the specified algorithm function. 
}

\usage{
gene2p(
  gwas,
  gts.source, 
  method = GATES, 
  cores = 1
)
SpD(ldmatrix, snps, p)
GATES(ldmatrix, snps, p)
}

\arguments{
  \item{gwas}{data.frame. Has to contain columns 'SNP', 'P' and either 'geneid' or 'genename' (when both are present, the one with the smaller index is used). A single SNP may occur with different genes, but always has to have the same p-value. }
  \item{gts.source}{vector(1). Can be a HapMap population identifier (numeric) to retrieve genotyes for, or an object of class \code{\link{snp.data}} holding genotypes, or a GenABEL (.gwaa) or LINKAGE / PLINK (.ped) genotype file, with existing corresponding .phe and .map files, respectively. See the \code{gts.source} argument in the \code{\link{getGenotypes}} function for the exact file format specifications. GenABEL format is fastest and recommended. The files might contain more SNPs than actually needed. The .map file has to be in --map3 format (columns CHR, SNP, BP). }
  \item{method}{function. Can currently be the SpD or GATES functions. See 'Details' for more information. }
  \item{cores}{integer. The number of parallel processes to use (cores = 1 uses no parallelization). }
  \item{ldmatrix}{numeric. A matrix of ld values, dimensions matching the \code{p} and \code{snps} arguments. }
  \item{snps}{character. A vector of SNP identifiers. }
  \item{p}{numeric. A vector of p-values matching the length of the \code{snps} argument. }
}

\details{
  The SpD and GATES functions calculate an aggregate p-value for a single gene/locus, based on multiple tested SNPs (e.g. an association test) that are assigned to that locus. To account for statistical dependence or independence between the multiple SNPs, linkage disequilibrium information 
  is required to decorrelate the tests. Further information about gene-based association tests and aggregate p-values and the specific GATEs and SpD methods can be found in the references section. 
  The gene2p function enables an application of the methods to multiple genes (or a genomewide dataset, with certain runtime limitations) and automatic calculation of the required LD amtrix for each gene.   
  Thus, conventionally the user will call the gene2p function to calculate a representative, gene-based p-value for a large set of genes that are annotated to multiple SNPs (e.g. by the \code{link{snp2gene}} function) as result of an association study.  
  Calculations may be time- and space demanding, so depending on the number of SNPs, it might be a good idea to divide the dataset by chromosome, use only intragenic SNPs or prune the the dataset (e.g. drop SNPs evenly).    
}

\value{
 data frame: The 'gwas' argument with a column 'gene.p' added. Contains the original number of rows. Column 'gene.p' may contain NA values when the gene-wise p-value could not be calculated for that row.
 As a side effect, the retrieved genotype data is deposited in files gene2p.gwaa, gene2p.phe, gene2p.ped and gene2p.map.   
}

\references{
  The GATES and SpD algorithms have been proposed by Miao-Xin Li et. al. in {http://dx.doi.org/doi:10.1016/j.ajhg.2011.01.019} and Dale Nyholt, respectively \url{http://dx.doi.org/doi:10.1086/383251}.
}

\examples{
  snps <- data.frame(SNP = c("rs188090", "rs172154", "rs759704"))
  snps$P <- runif(nrow(snps))^2
  \dontshow{
  ## offline usage: load and set buffer data for gene annotation
  setPostgwasBuffer(
    snps = read.table(system.file("extdata", "postgwas.buffer.snps", package = "postgwas"), header = TRUE), 
    genes = read.table(system.file("extdata", "postgwas.buffer.genes", package = "postgwas"), header = TRUE)  
  )
  }
  
  # offline LD annotation needs genotype files
  # hint: genotype data can also be preloaded using the getGenotypes function
  gwaafile <- system.file("extdata", "example.gwaa", package = "postgwas")
  
  gwas <- snp2gene.prox(snps, level = 0, use.buffer = TRUE)
  
  gene2p(gwas, gts.source = gwaafile)
  \dontshow{
    # cleanup
    clearPostgwasBuffer()
    rm(snps, gwaafile, gwas)
  }
}


