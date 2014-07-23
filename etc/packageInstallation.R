# package dependencies ====================================================
# setRepositories()
# setRepositories(ind=1:6)

depends  <- devtools::parse_deps(read.dcf("./DESCRIPTION")[, "Depends"])
suggests <- devtools::parse_deps(read.dcf("./DESCRIPTION")[, "Suggests"])
imports  <- devtools::parse_deps(read.dcf("./DESCRIPTION")[, "Imports"])

missing.depends.vanilla <- c('biomaRt', 'AnnotationDbi', 'org.Hs.eg.db', 'Rgraphviz', 'RBGL')
missing.suggests.vanilla <- c('topGO', 'VariantAnnotation', 'Rsamtools', 'org.Mm.eg.db', 'GO.db')

# if CRAN is the only mirror we will still be missing the following packages
# and use Bioconductor to install them
packages.CranOrInstalled <- c(available.packages()[,1], installed.packages()[, 1])
missing.depends <- depends[!depends[, 1] %in% packages.CranOrInstalled, 1]
missing.suggest <- suggests[!suggests[, 1] %in% packages.CranOrInstalled, 1]
missing.imports <- imports[!imports[, 1] %in% packages.CranOrInstalled, 1]

source("http://bioconductor.org/biocLite.R")
biocLite(missing.depends)
biocLite(missing.suggests)
biocLite(missing.imports)

bioconductorPackages <- c("biomaRt", "AnnotationDbi", "org.Hs.eg.db", "Rgraphviz", "RBGL")

# installing
devtools::install(".")






## procedure for README.md ============================================
# markdown: https://bitbucket.org/tutorials/markdowndemo

## Option 1 ####
setRepositories(ind = 1:6)
install.packages("postgwas")
library(postgwas)


## Option 2 ####
setRepositories(ind = 1:6)
install.packages("http://bitbucket.org/merns/postgwas/downloads/postgwas_1.11-2.zip", repos=NULL)
library(postgwas)


## Option 3 ####

# install devtools
install.packages(c("devtools", "rstudioapi"))
setRepositories(ind = 1:6)
devtools::install_bitbucket("postgwas", username="merns")
library(postgwas)
