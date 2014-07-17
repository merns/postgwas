% File man/postgwasBuffer.Rd
\name{postgwasBuffer}
\alias{postgwasBuffer}
\alias{getPostgwasBuffer}
\alias{setPostgwasBuffer}
\alias{clearPostgwasBuffer}

\title{Postgwas Buffer Variable Usage}

\description{
 This help page documents functions to retrieve and modify postgwas buffer data and contains explanatory text on the buffer data concept in postgwas.
}

\usage{
getPostgwasBuffer()
clearPostgwasBuffer()
setPostgwasBuffer(
    uselist = FALSE, 
    snps, 
    genes, 
    genes.regionalplot, 
    exons.regionalplot, 
    ld.regionalplot, 
    goterms
)
}

\arguments{
  \item{uselist}{list. This can be a list, where each list element corresponds to one buffer variable to be set. List elements . When this is not a list, the argument is ignored. The return value of \code{getPostgwasBuffer} is compatible with this argument, so buffer data can be quickly saved and restored with both methods. }
  \item{snps}{data frame. Sets to buffer variable for SNP annotation data to the corresponding data frame (execute the examples to see an example of such a data frame). When missing, the variable is unaltered. When NULL, the buffer variable is cleared (set to NULL). }
  \item{genes}{data frame. Sets to buffer variable for gene annotation data to the corresponding data frame (execute the examples to see an example of such a data frame). When missing, the variable is unaltered. When NULL, the buffer variable is cleared (set to NULL). }
  \item{genes.regionalplot}{data frame. Sets to buffer variable for gene annotation data in regionalplots to the corresponding data frame (execute the examples to see an example of such a data frame). When missing, the variable is unaltered. When NULL, the buffer variable is cleared (set to NULL). }
  \item{exons.regionalplot}{data frame. Sets to buffer variable for exon annotation data in regionalplots to the corresponding data frame (execute the examples to see an example of such a data frame). When missing, the variable is unaltered. When NULL, the buffer variable is cleared (set to NULL). }
  \item{ld.regionalplot}{list. Sets to buffer variable for ld annotation data in regionalplots to the corresponding list. List elements are named according to query SNP identifiers and contain a matrix of r square LD values for all SNPs in the region defined by the query SNP (see examples). When missing, the variable is unaltered. When NULL, the buffer variable is cleared (set to NULL). }
  \item{goterms}{data frame. Sets to buffer variable for GO term annotation data to the corresponding data frame (execute the examples to see an example of such a data frame). When missing, the variable is unaltered. When NULL, the buffer variable is cleared (set to NULL). }
}


\details{
Buffer data is stored in a environment named \code{postgwasBuffer} in the postgwas namespace (not accessible by the user with conventional means). 
The variables in that environment can be listed with \code{getPostgwasBuffer()}. They are returned as list and initially all set to NULL. 
Buffer data is only used when the argument \code{use.buffer} is set to TRUE in a function call (by default, is always FALSE).
When it is set to TRUE, downloaded annotation data is stored in the corresponding buffer variables as long as they are NUll, otherwise, existent data is used.
Buffer data can be set or modified manually using the \code{setPostgwasBuffer()} function. 
However, this is an advanced feature and only few checks are performed on the validity of the data. When setting \code{use.buffer = TRUE} always double-check whether the buffer data is valid for the specific function call. 
Which buffer variables are used by which function is listed in the according \code{use.buffer} argument documentation of the function. 
Buffer variables persist for a whole R session, but can be reset to NULL using the \code{clearPostgwasBuffer} function. 
}

\examples{
# initially, all buffers are empty
# there is one data frame for snps, genes, genes used only in regionalplots etc.
getPostgwasBuffer() 

\dontrun{
# downloaded data is stored in buffers when empty and the argument use.buffer is set
snp2gene.prox(
  snps = data.frame(SNP = c("rs172154", "rs759704"), 
  use.buffer = TRUE
)
head(getPostgwasBuffer()$snps)
}

# buffer data can also be set manually as a list containing all or some buffer variables
load(system.file("extdata", "bufferHS.RData", package = "postgwas"))
names(bufferHS)

# pass it to postgwas, overwriting previous data if existing
setPostgwasBuffer(bufferHS)

# show contents of all buffer data variables
lapply(getPostgwasBuffer(), head)

# also, single buffer elements can be set
# ld buffer for regionalplots has not been set yet 
getPostgwasBuffer()$ld.regionalplot

# load and set from file
load(system.file("extdata", "LDrs4846567region.RData", package = "postgwas"))
setPostgwasBuffer(ld.regionalplot = LDrs4846567)

# it is not a data frame but a list
# each list element is one ld matrix for each query SNP in the regionlplot call
names(LDrs4846567)
LDrs4846567[["rs4846567"]][1:6, 1:6]

# before starting a new analysis that relies on different annotation data, clear buffer 
clearPostgwasBuffer()

\dontshow{
  rm(bufferHS, LDrs4846567)
}
}