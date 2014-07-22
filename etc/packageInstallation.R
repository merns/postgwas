# package dependencies ====================================================
# setRepositories()
# setRepositories(1:6)

depends  <- devtools::parse_deps(read.dcf("./DESCRIPTION")[, "Depends"])
suggests <- devtools::parse_deps(read.dcf("./DESCRIPTION")[, "Suggests"])
imports  <- devtools::parse_deps(read.dcf("./DESCRIPTION")[, "Imports"])


bioconductorPackages <- c("biomaRt", "grid", "AnnotationDbi", "org.Hs.eg.db", "Rgraphviz", "RBGL")

# if CRAN is the only mirror we will still be missing the following packages
# and use Bioconductor to install them
missing.depends <- c('biomaRt', 'AnnotationDbi', 'org.Hs.eg.db', 'Rgraphviz', 'RBGL')
missing.suggests <- c('topGO', 'VariantAnnotation', 'Rsamtools', 'org.Mm.eg.db', 'GO.db')


source("http://bioconductor.org/biocLite.R")
biocLite(missing.depends)
biocLite(missing.suggests)


# installing
devtools::install_bitbucket("postgwas", username="merns", password="")
devtools::install(".")





## procedure for README.md ============================================
#
# markdown: https://bitbucket.org/tutorials/markdowndemo


# install devtools
install.packages(c("devtools", "rstudioapi"))

# install missing packages
biocPkgs <- c("biomaRt", "AnnotationDbi", "org.Hs.eg.db", "Rgraphviz", "RBGL")
source("http://bioconductor.org/biocLite.R")
biocLite(biocPkgs)

# install postgwas
devtools::install_bitbucket("postgwas", username="merns", password="")


## Alternative:
install.packages("https://bitbucket.org/merns/postgwas/downloads/postgwas_1.11-2.zip", repos=NULL)
