% File man/gwas2network.Rd
\name{gwas2network}
\alias{gwas2network}
\title{Highlight GWAS Results in Interaction Networks}

\description{
Visualize and analyze biological interactions between GWAS-associated genes. This helps to identify modules of sub-significantly associated genes in the GWAS that share a function, probably hinting on a polygenic origin of the investigated trait.
}

\usage{
gwas2network(
    gwas.mapped.genes, 
    network, 
    prune = "gwasonly",
    max.communities = 500, 
    vertexcolor.GO.overrep = "org.Hs.eg.db",
    vertexcolor.GO.regex = NULL,
    min.transparency = 0.25, 
    max.transparency = 0.75, 
    custom.layout = FALSE,
    edge.weight.fun = function(p1, p2, degree1, degree2, weight) {
      return(((-log10(p1) +1) * (-log10(p2) +1))^2 * weight)
    },
    edge.weight.collapse.fun = function(weights, labels) {
      return(mean(weights))
    },
    remove.superhubs = nrow(network) > 100,
    p.default = max(gwas.mapped.genes$P, na.rm = TRUE),
    png.maxedges = 2500, 
    file.verbosity = 2, 
    biomart.config = biomartConfigs$hsapiens, 
    use.buffer = FALSE, 
    cores = 1
)
}

\arguments{
  \item{gwas.mapped.genes}{data.frame. Has to contain columns "SNP", "P" and one or both of "genename" and "geneid" (one of these should match the identifiers in \code{network}. You might want to use the functions \code{\link{bm.id2name}} or \code{\link{bm.name2id}} for translation). When using non-GWAS derived data, the columns SNP and P can be set to dummy data (e.g. SNP = 1:n, P = 0.05) . Can optionally contain a "pheno" column to visualize serveral datasets in a single network (see Details). }
  \item{network}{data.frame or matrix. The first two columns each contain the ID or name of one of the interactors. Make sure that the identifiers used match one of the columns 'geneid' or 'genename' in gwas.mapped.genes (automatically selects one of both columns based on which has the higher overlap with identifiers in the network argument). Can also contain a "label" column to display edge labels and a weight column (see details and the \code{edge.weight.fun} argument). Loops and duplicates are automatically removed from the network. See also \code{\link{getInteractions.path}} function. }
  \item{prune}{character. Can be NULL, "gwasonly", "connected" and "shared". When not NULL, the network will be reduced: "gwasonly" keeps only edges between genes from the \code{gwas.mapped.genes} argument, "shared" keeps all between \code{gwas.mapped.genes} and additionaly genes that form edges with two distinct genes from \code{gwas.mapped.genes}, and connected" keeps additionally genes that have an edge with any one gene from the \code{gwas.mapped.genes}.}
  \item{max.communities}{integer(1). When > 1, runs \code{\link[igraph]{spinglass.community}} on all connected components and stores the resulting communities in separate plots. Creates at most max.communities communities per connected component (by setting the 'spins' parameter in \code{\link{spinglass.community}}). Maximum value is 500. When vertex colorization by GO term overrepresentation analysis is selected, each community will be tested for overrepresentation tested separately. When <= 1, plots the connected components in separate plots, reusing the vertex colorization of the whole graph. }
  \item{vertexcolor.GO.overrep}{character(1). Has to be the package name of an installed GO annotation package (e.g. org.Hs.eg.db, see also \url{http://www.bioconductor.org/packages/2.10/data/annotation/}) When this is set, an overrepresentation analysis (fisher exact test and 'weight01' decorrelation algorithm from the topGO package) of network genes in gene ontology terms will be performed on the whole graph, and vertices colorized corresponding to the three most significant terms. Further, a legend with p-value and color code for the terms is added. Communities will get a separate overrepresentation analysis for each community. Connected components inherit the vertex color scheme from the single whole-graph overrepresentation analysis. \emph{Is only considered when the vertexcolor.GO.regex argument is NULL}. When both are NULL, vertex colors are black.}
  \item{vertexcolor.GO.regex}{list. A named list of regular expressions. The names have to be valid color names, the values are parsed as regular expressions. Gene ontology terms will be downloaded for each gene, and genes colorized that match their annotated terms with an expression. Genes that match multiple expressions are colored brown. Has precedence over the vertexcolor.GO.overrep argument. When both are NULL, vertex colors are black.}
  \item{min.transparency}{numeric, in [0,1]. Sets the minimum transparency of vertices. Also affects edges at min.transparency^1.8. Vertex p-values define the transparency degree, and this argument defines the lower bounds for the transparency scale. }
  \item{max.transparency}{{numeric, in [0,1]. Sets the maximum transparency of vertices. Also affects edges at max.transparency^1.8. Vertex p-values define the transparency degree, and this argument defines the upper bounds for the transparency scale. }}
  \item{remove.superhubs}{boolean. When TRUE, all vertices are removed that are not in gwas.mapped.genes and satisfy degree(vertex) > mean(degree(network))^2. }
  \item{custom.layout}{boolean(1). When TRUE, opens a window with draggable vertices to create a custom layout of the plotted graph. Otherwise, the Fruchtermann-Reingold algorithm is used which usually produces very good results. Note that the vertex colorization, label to vertex positions and vertex size are not preserved in the manual layout window.}
  \item{p.default}{numeric(1). Has to be in [0,1]. When the network contains vertices that are not present in gwas.mapped.genes, they will get this default p-value assigned. This affects the displayed size of these vertices and the edge weights during community detection. By default, this is the largest p-value (or gene-based p-value) that is listed in \code{gwas.mapped.genes}. }
  \item{edge.weight.fun}{function. Calculates the weight of edges based on vertex weight, degree and preexisting edge weights. This is a parameter needed for community detection, thus only used when detect.comunities = TRUE. See details.)}
  \item{edge.weight.collapse.fun}{function. When there is a preexisting \code{weight} column given in the network argument, and there are duplicate / multiple edges, this function is applied to calculate a combined weight (is applied before the inclusion of vertex P, i.e. edge.weight.fun). The arguments of this function get passed the vector of weights of the duplicate edges, and the vector of labels (labels is NULL when no labels are present in the \code{network}argument).}
  \item{png.maxedges}{integer(1). Only graphs with edgecount <= png.maxedges will be plotted to a png file.}
  \item{file.verbosity}{integer. Can be 0 for no output to files, 1 for image and pdf files of graphs, 2 for additional graph text files (multi annotation data, graph data, global parameter settings etc) and 3 for additional overrepresentation data.}
  \item{biomart.config}{A biomart configuration list. See \code{\link{biomartConfigs}}}
  \item{use.buffer}{boolean. When TRUE, buffers the data downloaded for gene name to ID mapping in the \code{genes} buffer variable. See \code{\link{postgwasBuffer}} for more information on buffer variables. When vertegxcolor.GO.regex is set, GO term annotations for each gene is stored in the \code{goterms} buffer variable. When this variables are already set, data is not downloaded but used from there instead. This facilitates the possibility to re-use data, alter data or provide custom data. When setting use.buffer = TRUE, make sure that pre-existing buffer data is current!}
  \item{cores}{integer. Specifies the allowed number of parallel processes. Requires the package 'multicore' to be installed when more than 1 process should be used for computation. Consumes a multiplicity of memory according to number parallel processes started. WORKS CURRENTLY ONLY ON LINUX SYSTEMS. Default is 1 (multithreading deactivated)}
}


\details{
This function displays the given network and highlights all genes from the argument gwas.mapped.genes and their connecting edges using transparency effects. 
Beside returning the resulting graph object, several files will be generated containing diverse output data: 
\describe{
  \item{gwas2network.pdf}{Plot of the complete network. Depending on the number of edges, an accompanying png file is produced as well.}
  \item{gwas2networkCommunities*.pdf}{By default, a graph partitioning algorithm is applied to additionally decompose connected components of the network into (functional) modules (see below for details). These are afterwards classified by GO term overrepresentation analysis and labelled appropriately. Additional png output is as above.}
  \item{gwas2networkConnectedComponents*.pdf}{When graph partitioning is turned off (by setting max.communities <= 1), only connected components are plotted separately. Additional png output is as above.}
  \item{gwas2networkParams.txt}{Contains most of the function argument values.}
  \item{gwas2networkMultiAnnotations.txt}{A list of ambigouos SNP to gene annotations (multiple genes that inherit the p-value from a single association). These are displayed as vertices with a white cross in the plot (see below). When multiple phenotypes are present, one file will be constructed for each phenotype.}
  \item{gwas2networkGraphEdges.csv}{The edge data of the graph after processing by gwas2network in a tab delimited format.}
  \item{gwas2networkGraphVertices.csv}{Vertex data.}
  \item{gwas2networkGOoverrep.csv}{The complete list of significant terms (uncorrected) of the GO overrepresentation analysis in the whole graph.}
}

Data for the 'network' argument can be retrieved using the separate \code{getInteractions} functions of the package. The network is normalized and pruned according to the \code{prune} argument. Normalization includes collapsing of multi-edges to unique edges (concatenates labels when edge lables are used), removal of loops and occasionally removal of hubs with overproportional high number of connections. The latter only applies to non-GWAS genes and makes sense because promiscuously interacting non-GWAS genes usually do not add much specifity to a module of interacting GWAS genes, which is what we would like to find in the network. For small networks, superhub removal should be disabled (default), because vertices with a high degree in a small network do not need to be hubs in the entire interaction universe and should thus be retained. Negative weights to edges are principally allowed, but not recommended. For negative weights the clustering argument will be inverted (i.e. the algorithm tries to place negative edges between communities, see \code{link[igraph]{spinglass.community}})

Vertex data: The "P" column in gwas.mapped.genes determines the size of the vertices (the -logp10 id calculated ad the larges of these values set to a size of 1). May not be negative. 
Usage of names vs IDs: When IDs are used as base identifier, all analyses will operate on IDs until visualization. Then, IDs are translated to names (using the biomart config) to ensure readability of the graph. When names are used as base identifier, no translation takes place except for the GO overrepresentation analysis, where entrez gene ids have to be used. 
When a "pheno" column is present, a unique vertex shape is assigned to each phenotype (currently capped to 3 different phenotype categories) and displayed accordingly. When a vertex (gene) falls into several phenotype categories, the vertex name is displayed boldface and italic, and all vertex shapes for all phenotypes are plotted on top of each other. The best P-value vertex is plotted first with the label and graph edges assigned to it, and the other vertices are plotted on top of that, but the shape is reduced to its frame and having a slightly differing color intensity. This way it is possible to identify the p-value of the multi-phenotype vertices.

There can also be the case that a SNP is annotated to more than one gene. These genes will have a white cross-mark plotted on their vertices, indicating the uncertainity of annotation. This is independent of the phenotypes, i.e. ambiguous annotation is detected only within each phenotype subset. This feature can be 'deactivated' by assigning unique dummy SNP identifiers in the gwas.mapped.genes argument (the SNP column is used only for ambiguous annotation detection).

The graph partitioning method (\code{link{igraph::spinglass.community}}) considers edge weights to find densely connected modules (communities). To obtain modules that reflect a dysregulated molecular pathway through multiple associated GWAS genes, we should define edge weights that are higher for a pair of trait-associated genes than unassociated genes. The exact edge weight function can be user specified. The default weight function is designed with the following considerations: We use the log10 of the vertex p-values divided by a penalty for highly connected vertices. The penalty decreases logarithmically for large degree (e.g. does not change much for 100 and 200 connections) and involes addition of a constant term to decrease the impact of penalty for few connections. This makes sense because formation of a module should comprise genes that interact specifically, which normally decreases for genes with more than ~20 connections. The base of the logarithm and the constant shift are derived empirically and might need slight change for different networks (e.g. you might want to increase the shift when the network is generally densely connected or not hub-based). A custom edge weight via a \code{weight} column in the \code{network} parameter is also considered in the function and is used as a multiplicator for the total edge weight by default. When the \code{weight} column does not exist it is set to 1, and so are NA values in that column. Duplicate edges will be collapsed to the mean weight of all duplicate edges.

To emphasize direct connections between genes associated in the GWAS, the function return value could be squared,
but in general this is not recommended (tends to break modules apart when they use intermediate vertices)
}

\value{
  Returns a graph as igraph object. 
  Plots further the resulting graph and its communities to file, along with several text summary files (see 'Details').
}

\seealso{
  \code{\link{biomartConfigs}}, \code{\link{snp2gene}}, \code{\link{igraph}}, \code{\link{getInteractions.proteincomplex}}, \code{\link{getInteractions.domains}}, \code{\link{getInteractions.GO}}
}

\examples{

\dontshow{
  ## preload packages
  require(igraph, quietly = TRUE)
  
  ## offline usage: load and set buffer data
  setPostgwasBuffer(
    snps = read.table(system.file("extdata", "postgwas.buffer.snps", package = "postgwas"), header = TRUE), 
    genes = read.table(system.file("extdata", "postgwas.buffer.genes", package = "postgwas"), header = TRUE)  
  )
}

# read example gwas dataset
gwas.result <- read.table(
                  system.file("extdata", "example.gwas1", package = "postgwas"), 
                  header = TRUE
               )
snps <- gwas.result[gwas.result$P < 1e-03, ]

# we annotate all neighbor genes for simplicity and construct a network on these
gwasgenes.prox <- snp2gene.prox(snps, use.buffer = TRUE)


# get network data
\dontrun{
  # download and use combined protein interaction and domain interaction data
  network.ppi <- getInteractions.proteincomplex(filter.ids = gwasgenes.prox$geneid)
  network.ppi$label <- "ppi"
  network.dom <- getInteractions.domains(
                   filter.ids = gwasgenes.prox$genename, 
                   max.occurence = NULL
                 )
  network.dom$label <- network.dom$domain.name
  network <- rbind(network.ppi, network.dom[, c("geneid.x", "geneid.y", "label")])
}

\dontshow{
network <- read.table(system.file("extdata", "postgwas.interaction.download", package = "postgwas"), header = TRUE)
}

# run gwas2network without vertex color
# on windows the on screen device does not produce good results, 
# set file.verbosity > 0 to plots graph + communities cleanly to files
# note that entrez IDs are passed as network argument, but names are plotted
net <- gwas2network(
  gwasgenes.prox,
  network = network,
  prune = "shared",
  vertexcolor.GO.overrep = NULL, 
# uncomment this line to colorize by regex matching GO terms  
#  vertexcolor.GO.regex = list(red = "binding", blue = "membrane"),
  file.verbosity = 0, 
  use.buffer = TRUE
)
# display the returned graph
gwas2network.plot(net, device = options("device")$device)

\donttest{

  ###### comparative network between two studies plus overrepresentation analysis ######

  # construct random association data of a second study 
  # (uncomment the following line to re-generate random IDs)
  # random.geneids <- as.numeric(
  #   mappedkeys(org.Hs.egGENENAME)[floor(runif(100, 1, length(as.list(org.Hs.egSYMBOL))))])
  random.geneids <- c(477, 4751, 9899, 10316, 10776, 29901, 80124, 
                      221395, 93973, 157313, 5536, 4735, 55089)
  gwasgenes.random <- data.frame(
    # dummy SNP names, no multi annotations
    SNP = 1:length(random.geneids), 
    geneid = random.geneids,
    P = c(runif(length(random.geneids) -1, 0.00000001, 0.02)^2, 2*10^-6)
  )

  gwasgenes.prox$pheno <- "study1"
  gwasgenes.random$pheno <- "study2"
  gwasgenes.combined <- rbind(gwasgenes.prox[, c(6,1,9,12)], gwasgenes.random)

  net.GO <- getInteractions.GO(as.vector(gwasgenes.combined$geneid), toFile = NULL)
  net.GO <- net.GO[net.GO$weight > 0.75, ]

  # this produces several result files including plots in the current directory!
  net <- gwas2network(
    gwasgenes.combined,
    network = net.GO,
    file.verbosity = 0
  )
  gwas2network.plot(net, device = options("device")$device)
  
  # when file.verbosity >= 2, we can restore graph data from text files to do further study
  # v <- read.table("gwas2networkGraphVertices.csv", header = TRUE)
  # e <- read.table("gwas2networkGraphEdges.csv", header = TRUE)
  # g <- graph.data.frame(e, vertices = v, directed = FALSE)

  \dontshow{
    # cleanup
    rm(gwasgenes.combined, gwasgenes.random, net.GO, random.geneids)
  }
}
\dontrun{

  ###### example for yeast ######
  gwas.result <- read.table(
                   system.file("extdata", "example.gwas3.xz", package = "postgwas"), 
                   header = TRUE
                 )
  snps <- gwas.result[gwas.result$P < 1e-03, ]
  gwasgenes.prox <- snp2gene.prox(snps, biomart.config = biomartConfigs$scerevisiae)
  
  network <- getInteractions.path(
               filter.ids = gwasgenes.prox$geneid, 
               biomart.config = biomartConfigs$scerevisiae
             )
  
  # use custom GO colors, because GOSim is by default configured for human data only
  net <- gwas2network(
    gwasgenes.prox,
    network = network,
    vertexcolor.GO.regex = list(red = "transcription", blue = "membrane"),
    biomart.config = biomartConfigs$scerevisiae
  )
  gwas2network.plot(net, device = options("device")$device)
  
  
  # now if we want to do it with GO overrepresentation analysis, 
  # we have to reconfigure GOSim and the IDs used...
  
  # configure GOSim for yeast
  # install yeast GO data from bioconductor 
  # (http://www.bioconductor.org/packages/2.10/data/annotation/)
  if(is.na(match("org.Sc.sgd.db", installed.packages()[, "Package"]))) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("org.Sc.sgd.db")
  }
  
  # GOSim does not take entrez IDs for yeast 
  # (this is because of the altered configuration in the org.Sc.sgd.db package...)
  # we can overcome this by changing the biomart config to use the appropriate IDs
  myconfig <- biomartConfigs$scerevisiae
  myconfig$gene$attr$id <- "external_gene_id"
  gwasgenes.prox <- snp2gene.prox(snps, biomart.config = myconfig)
  
  we can now also use a GO network if desired
  library(org.Sc.sgd.db)
  network <- getInteractions.GO(
              gwasgenes.prox$geneid, 
              GOpackagename = "org.Sc.sgd.db", 
              toFile = NULL
            )
  
  net <- gwas2network(
    gwasgenes.prox,
    network = network, 
    vertexcolor.GO.overrep = "org.Sc.sgd.db", 
    max.communities = -1,
    biomart.config = myconfig
  )
  gwas2network.plot(net, device = options("device")$device)
}

\dontshow{
  clearPostgwasBuffer()
  rm(gwasgenes.prox, gwas.result, net, network, snps)
}

}
