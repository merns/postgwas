% File man/gwas2network.plot.Rd
\name{gwas2network.plot}
\alias{gwas2network.plot}
\title{Automatically Scaled Graph Plot}

\description{
This is a replacement for the default plot function of igraph. It mainly offers simplified scaling properties, scaling the size of the medium proportional to the graph size and maintaining all graph elements at proper size. 
Compared to the default igraph plot function, it is truncated in other generic functionality and a bit slower. 
}

\usage{
gwas2network.plot(
  g = g,
  layout = layout.fruchterman.reingold(g, area = vcount(g)^3),
  legend.autolayout = TRUE, 
  vertex.to.edge.space = 1.3, 
  device = pdf, 
  ...
) 
}

\arguments{
  \item{g}{An igraph graph object. Requires vertex attibutes 'weight', 'color' and 'label' and the edge attribute 'color'. Can contain additional attributes (see Details section).}
  \item{layout}{matrix. Vertex x and y positions are listed here, refer to the igraph layout functions for detailed information, e.g. the 'Value' section of \link[igraph]{layout.fruchterman.reingold}.}
  \item{legend.autolayout}{boolean. When TRUE, creates a separate area and assigns all legend vertices to it. Legend vertices are defined as having the (sub)string 'postgwas' in their name. When FALSE, plots all vertices to the original coordinates given by the 'layout' argument.}
  \item{vertex.to.edge.space}{numeric. The space between vertices and their connecting edges, it is normally visually appealing to have a small gap between vertex and edge.}
  \item{device}{function. A device function that directs the plot data to a device, e.g. \code{pdf}, \code{X11}, ... }
  \item{...}{These additional parameters are passed to the \code{device} function. Can for example be \code{file} to specify a filename, see the appropriate device function documentation for valid arguments. Can also specify an argument 'v.weight.max' for vertex size normalization (the largest size will then correspond to v.weight.max instead of the maximum weight in the given graph).}
}

\details{
Compared to the default igraph plot, this function offers improved label positioning, cute edge endings for non-directed graphs and automatic scaling of all graph elements with respect to the graphs size. 
Its main use is to plot specific biological interaction networks and requires special graph attributes, and offers little customization properties. Draws always undirected networks.
The vertex attribute weight has to be a positive number and determines the size of the vertices. 
There can be extra attributes:
Edges can have a label attribute that is drawn at the middle of each edge.
Vertices can have the following optional attributes:
marked: boolean vertex attribute that plots a white cross on vertices where marked == TRUE
categories: vertices can optionally fall in categories. Currently three categories are allowed which correspond to vertex shapes. These are 'circle', 'square' and 'diamond'. There is one main category. Multiple categories per vertex are allowed. In this case, multiple vertices will be plotted on top of each other (the edges are only connected to the main category vertex). Several vertex attributes are used to determine main category and multiple categories. The main category is listed in a vertex attribute 'category.main' that contains as value the allowed category string 'circle', 'square' or 'diamond' (unknown values are set to one of these automatically). For each category that is used, there has to be an additional vertex attribute which is named as the category itself and is boolean. These additional attributes are mandatory even when there are no multi-categories (all FALSE). When a vertex falls into several categories, such a column value can be set to TRUE to denote multiple categories. Lastly, we have for each category the weight and marked status given in the vertex attribues circle.weight, square.weight and diamond.weight and accordingly circle.marked, square.marked and diamond.marked (also mandatory). These are considered for plotting the multi-category with the according parameters. Labels of multi-category vertices are printed in bold-italic. 
Size of different vertex shapes is always area-proportional with the weight argument (not same side length / radius).
}

\value{
 Is invoked for its side effect of plotting a graph to a specific device. 
}

\seealso{
  \code{\link{gwas2network}}, \code{\link{device}}
}
