% File man/getInteractions.Rd
\name{getInteractions}
\alias{getInteractions.path}
\alias{getInteractions.proteincomplex}
\alias{getInteractions.domains}
\alias{getInteractions.GO}
\title{Retrieve Biological Interaction Network Data}

\description{
 A collection of functions to download or calculate network data for biological interaction (or relatedness, co-occurence etc). The getInteractions.GO function uses functionality from the GOSim source code by Holger Froehlich (see references). Parts of that code have been included and adapted for parallelization.  
}

\usage{
getInteractions.path(
    filter.ids = "", 
    biomart.config = biomartConfigs$hsapiens,
    ds = bm.init.path(biomart.config), 
    toFile = "postgwas.interaction.download"
)
getInteractions.proteincomplex(
    filter.ids = "", 
    biomart.config = biomartConfigs$hsapiens,
    ds = bm.init.proteincomplex(biomart.config), 
    toFile = "postgwas.interaction.download"
)
getInteractions.domains(
  filter.ids = "", 
  biomart.config = biomartConfigs$hsapiens,
  filter.type = biomartConfigs$hsapiens$gene$filter$name,
  ds = bm.init.genes(biomart.config),
  min.occurence = NULL, 
  max.occurence = if(length(filter.ids) < 20)  
                    NULL 
                  else 
                    length(filter.ids) / (logb(length(filter.ids), 3) +1),
  toFile = "postgwas.interaction.download"
)
getInteractions.GO(
  filter.ids, 
  GOpackagename = "org.Hs.eg.db", 
  ontology = "BP", 
  toFile = "postgwas.interaction.download", 
  cores = 1, 
  ...
)
}

\arguments{
  \item{filter.ids}{vector. A list of identifiers that that are used to retrieve network data. The type of IDs is specified by additional parameters (can also be symbols). See details.}
  \item{biomart.config}{list. Specifies field names for biomart data retrieval. See \code{\link{biomartConfigs}}.}
  \item{filter.type}{character. A biomart filter field name, that specifies the type of IDs listed in the argument \code{filter.ids}.}
  \item{ds}{object. Biomart dataset to use. Is initialized automatically, but can also be explicitely passed to reduce runtime and network load in multiple calls.}
  \item{min.occurence}{numeric(1). Specifies the minimum total number of occurences for a domain to be considered for interaction. This can occasionally be useful to reduce the specifity of modules and complexity of the network., i.e. only domains with common functions will form modules. When NULL, this is deactivated (default). }
  \item{max.occurence}{numeric(1). Specifies the maximal total number of occurences for a domain to be considered for interaction. This can be used to avoid very common domains to create giant modules while increasing the functional specifity of modules. When NULL, this is deactivated. Is by default a log fraction of the number of vertices.}
  \item{GOpackagename }{character(1). Has to be the package name of an installed GO annotation package (e.g. org.Hs.eg.db, see also \url{http://www.bioconductor.org/packages/2.10/data/annotation/}). GO term information from this package will be used for similarity calculation. }
  \item{ontology}{character(1). Can be "BP", "MF" and "CC". Is passed to the \code{setOntology} function in GOSim, see the corresponding documentation.}
  \item{toFile}{character(1). A file/pathname where the downloaded interaction data is saved in tab-delimited format. When NULL, does not write a resultfile.}
  \item{cores}{integer. The number of parallel processes to use (cores = 1 uses no parallelization). }
  \item{...}{Optional arguments that are passed to the \code{getGeneSim} function of the GOSim package. Could be in particular the \code{similarity} argument from that function. }
}

\details{
\describe{
  \item{getInteractions.path}{Uses biomart (REACTOME) to determine pathway annotations for each gene. Common pathways between genes form interactions. Uses Entrez gene IDs by default.}
  \item{getInteractions.proteincomplex}{Uses biomart (REACTOME) to determine protein interactions (from known complexes) between query genes. Uses Entrez gene IDs by default.}
  \item{getInteractions.domains}{Uses biomart to annotate (superfamily) domains to each gene. The edges are formed by common domains between proteins. Uses gene symbols by default.}
  \item{getInteractions.GO}{Calculates relations between the query genes based on similarity of their GO term architecture (using the GOSim package). This is time-intensive and takes seveal hours for ~1000 genes. Note that interactions at arbitrary strength are returned (consider 'weight' column), these should manually be filtered for a minimum similarity of e.g. weight > 0.75. GO term retrieval depends on bioconductor annotation packages. The human (org.Hs.eg.db) package is installed by default. When switching to a different organism by installing the appropriate package (\url{http://www.bioconductor.org/packages/2.10/data/annotation/}), the type of IDs used may change. The IDs used are then specified by the annotation package. See the example section of \link{gwas2network} how to overcome such an ID problem exemplified with yeast.} 
}  
}

\value{
  A data frame with columns
  \describe{
    \item{getInteractions.path:}{"geneid.x", "geneid.y" and "label". The first two columns contain Entrez gene IDs by default, and the third the shared pathway names. }
    \item{getInteractions.proteincomplex:}{"geneid.x", "geneid.y" and "label". The first two columns contain Entrez gene IDs by default, and the third the shared protein complex names. }
    \item{getInteractions.GO:}{"geneid.x", "geneid.y" and "weight". The first two columns contain gene IDs / names according to the annotation package configuration, currently entrez gene IDs for the human package org.Hs.eg.db. The weight column contains the similarity measure between the genes calculated by GOSim. }
    \item{getInteractions.domains:}{returns columns "genename.x", "genename.y", "geneid.x", "geneid.y", "domain.name", the latter specifying the shared domain that defines the edge, and IDs and names according to the biomartConfig used.}
  }
}

\seealso{
  \code{\link{biomartConfigs}}
}

\references{
  The GOSim source code used by the getInteractions.GO function is described in the following publication: Holger Froehlich, N. Speer, A. Poustka, Tim Beissbarth (2007). GOSim - An R-Package for Computation of Information Theoretic GO Similarities Between Terms and Gene Products. BMC Bioinformatics, 8:166

  The Reactome database used by the getInteractions.path function can be accessed at \url{http://www.reactome.org}
}
