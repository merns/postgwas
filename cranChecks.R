
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

# devtools::check(cran=TRUE, cleanup=FALSE)
# uses roxygen2!
# dont use devtools::check, as it overwrites the NAMESPACE file




## Manual CRAN checks ##############

# .Rbuildignore adjustment (Remark: Remove double backslashes from output of glob2rx())
# glob2rx("postgwas_1.11-1.tar.gz")

cranBuild <- paste0("\"", Sys.getenv("R_HOME"), "/bin/R.exe\" CMD build .")
shell(cranBuild)
cranCheck <- paste0("\"", Sys.getenv("R_HOME"), "/bin/R.exe\" CMD check --as-cran postgwas_1.11-1.tar.gz")
shell(cranCheck)

vignette(tables)

