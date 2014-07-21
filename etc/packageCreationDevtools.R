# Package creation with "devtools" ===============================================

# create vignette files in "./inst/doc" for instance for when vignette building is
# disabled in DESCRIPTION with "BuildVignettes: FALSE"
devtools::build_vignettes()   # takes long time

# builds (bundles) package to *.tar.gz and puts it in parent directory
devtools::build(vignettes=FALSE)  #calls R CMD build


# builds and checks the package
# usually uses roxygen2, which would break our non-roxygen2 documentation
devtools::check(cran=TRUE, document=FALSE, doc_clean=FALSE)

# producing a binary packages sometimes fails for no known reason when testing to load i386 libs
# use manual version in "packageCreationCRAN.R" to have a save way producing binary packages or use build_win() function.
devtools::build(binary=TRUE, vignettes=FALSE) #calls R CMD install

# build binary with win-builder.r-project.org
# keep in mind that build is send to current maintainer in DESCRIPTION file
devtools::build_win()

# release package to CRAN with the help of devtools function
devtools::release()
