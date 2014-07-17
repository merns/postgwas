% File man/ytracks.regionalplot.Rd
\name{ytracks.regionalplot}
\alias{ytracks.regionalplot}
\title{Generate a Base Tracks Object for Regionalplots}

\description{
  This is a function to construct a data frame listing y axis bounds for tracks in the regionalplot function. 
}

\usage{
  ytracks.regionalplot(
    ylim.upper, 
    ylim.space, 
    plot.genes = TRUE, 
    plot.ld = TRUE, 
    plot.variants = FALSE, 
    var.options
  )
}

\arguments{
  \item{ylim.upper}{numeric(1). The height of the plot starting at 0, i.e. the maximum scale of the p-value graph.}
  \item{ylim.space}{numeric(1). A default value for gaps between tracks.}
  \item{plot.genes}{boolean(1). Whether to add information for a gene track.}
  \item{plot.ld}{boolean(1). Whether to add information for an LD track.}
  \item{plot.variants}{boolean(1). Whether to add information for a rare variant track.}
  \item{var.options}{list. Only needed when plot.variants = TRUE, should then be the same list as passed to the \code{\link{regionalplot}} function.}
}

\details{
  See also the last example of the \code{\link{regionalplot}} function, where the use of the tracks function is outlined.
}

\value{
  A data frame with each track forming a row and the columns listing the properties (see examples).
}

\seealso{
  \code{\link{regionalplot}}
}

\examples{
  ytracks.regionalplot(ylim.upper=12, ylim.space=0.1)
}