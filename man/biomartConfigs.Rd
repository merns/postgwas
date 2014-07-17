% File man/biomartConfigs.Rd
\name{biomartConfigs}
\alias{biomartConfigs}
\title{Configuration Lists for Biomart Access}

\description{
  This is a list. For each organism, it contains sublists with appropriate field names for biomart queries. 
}

\details{
Several postgwas functions rely on SNP, gene and exon data retrieval from biomart. All these functions use one type of configuration list that contains the necessary biomart field names. Customizations and adaptation to changes in the field names of a biomart database can thus be achieved by altering those configuration list. 

Structure: The first level contains configuration lists for diverse organisms. The available species can be listed by stating names(biomartConfigs). 
On the second level, it contains the configuration itself. i.e. the necessary field names.  The structure can be investigated by stating \code{names(unlist(biomartConfigs$hsapiens))}. The values of these elements can be changed as long as they are valid in the corresponding biomart database. For example, setting \code{biomartConfigs$gene$attr$name <- "uniprot_swissprot_accession"} will plot or assign accession numbers instead of gene symbols. The \code{\link{listAttributes}} function from the biomart package can be used to display valid field names. For some organisms, exon information is not available by definition. When one of the 'exon' attribute elements is set to NA, exon information will be ignored. 
}

\examples{
  
  # extract a configuration list
  myconfig <- biomartConfigs$drerio
  
  # alter the IDs used - first find which are available
  \dontrun{
    attr.avail <- listAttributes(useMart(
                    biomart = myconfig$gene$mart, 
                    dataset = myconfig$gene$dataset
                  ))
    # find the biomart field name for protein accession numbers
    attr.avail[grep("accession", attr.avail$description, ignore.case = TRUE), ]
  }
  
  # customize the configuration list to return accession numbers
  myconfig$gene$attr$name <- "uniprot_swissprot_accession"
  
  # also, use ENSEMBL IDs for filtering and retrieval
  myconfig$gene$attr$id <- "ensembl_gene_id"
  myconfig$gene$filter$id <- "ensembl_gene_id"
  
  # show elements of the new configuration list (genes only)
  myconfig$gene
  
  # cleanup
  rm(myconfig)
  
}
