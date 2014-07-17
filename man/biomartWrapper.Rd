% File man/biomartWrapper.Rd
\name{bm.id2name}
\alias{bm.name2id}
\alias{bm.id2name}
\title{Conversion Between Gene IDs And Names}

\description{
 Uses biomart to convert from gene names to IDs and back.  
}

\usage{
bm.id2name(
  filter.ids = NULL,
  config = biomartConfigs$hsapiens, 
  ds = bm.init.genes(config), 
  use.buffer = FALSE
)
bm.name2id(
  filter.names = NULL,
  config = biomartConfigs$hsapiens, 
  ds = bm.init.genes(config), 
  use.buffer = FALSE
)
}


\arguments{
  \item{filter.ids}{Either a vector containing IDs, or a data frame with an appropriate colum. The type of IDs should match the biomart config field gene$filter$id.}
  \item{filter.names}{Either a vector containing gene names, or a data frame with an appropriate colum. The type of namess should match the biomart config field gene$filter$name.}
  \item{config}{A list containing values that are needed for biomaRt connection and data retrieval. See \code{\link{biomartConfigs}}.}
  \item{ds}{Optionally, a biomart connection object for re-use scenarios to prevent connection overhead.}
  \item{use.buffer}{Uses the buffer variable \code{genes} for annotation when it exists and contains valid data in biomart format. See \code{\link{postgwasBuffer}} for more information on buffer variables. Does not write that variable.}
}

\details{
The filter data should not contain NAs.

# all names are unique
# only uses but does not set the buffer variable when use.buffer = TRUE
# use.buffer causes a merge with buffer data on all common columns

}

\value{
  Returns a data frame with an additional column that is named as the biomart config field for Ids (gene$attr$id) or names (gene$attr$name), respectively.
}

\seealso{
  \code{\link{biomartConfigs}}
}
