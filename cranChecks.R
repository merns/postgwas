## update devtools ##############

devtools::build_github_devtools()
#### Restart R before continuing ####
install.packages("./devtools.zip", repos = NULL)
# Remove the package after installation
unlink("./devtools.zip")



## package dependencies ##############

# depends
install.packages(c('data.table', 'biomaRt', 'GenABEL', 'igraph', 'gtools', 'gdata', 'coin', 'AnnotationDbi', 'org.Hs.eg.db', 'Rgraphviz', 'RBGL'))
# suggests
install.packages(c('topGO', 'VariantAnnotation', 'Rsamtools', 'org.Mm.eg.db', 'xtable', 'genetics', 'GO.db'))

missing.depends <- c('biomaRt', 'AnnotationDbi', 'org.Hs.eg.db', 'Rgraphviz', 'RBGL')
missing.suggests <- c('topGO', 'VariantAnnotation', 'Rsamtools', 'org.Mm.eg.db', 'GO.db')
missing.namespace <- c("BatchJobs")
install.packages(missing.namespace)


# setRepositories()
# setRepositories(1:6)

source("http://bioconductor.org/biocLite.R")
biocLite(missing.depends)
biocLite(missing.suggests)




## CRAN checks ##############

devtools::build(vignettes=FALSE) # source
#devtools::build_vignettes()   # takes long time

# builds and checks the package
# usually uses roxygen2, which would break our non-roxygen2 documentation
devtools::check(cran=TRUE, document=FALSE, doc_clean=FALSE)

# build binary
devtools::build(binary=TRUE)
# build binary with win-builder.r-project.org
devtools::build_win()


devtools::install_bitbucket("postgwas", username="merns", password="")
vignette("postgwas")


## Manual CRAN checks ##############

# .Rbuildignore adjustment (Remark: Remove double backslashes from output of glob2rx())
# glob2rx("postgwas_1.11-1.tar.gz")

setwd(paste0(getwd(), "/.."))

# build .tar.gz
cranBuild <- paste0("\"", Sys.getenv("R_HOME"), "/bin/R.exe\" CMD build postgwas")
shell(cranBuild)
# check .tar.gz (otherwise source tree would be tested)
cranCheck <- paste0("\"", Sys.getenv("R_HOME"), "/bin/R.exe\" CMD check --as-cran postgwas_1.11-1.tar.gz")
shell(cranCheck)
# usual way to build and install package
cranInstall <- paste0("\"", Sys.getenv("R_HOME"), "/bin/R.exe\" CMD INSTALL ../postgwas_1.11-1.tar.gz")
shell(cranInstall)
cranInstall <- paste0("\"", Sys.getenv("R_HOME"), "/bin/R.exe\" CMD INSTALL --build ../postgwas_1.11-1.tar.gz")
shell(cranInstall)

setwd(paste0(getwd(), "/postgwas"))




