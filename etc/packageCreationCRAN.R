## Package creation according to CRAN ==============================================
#
# http://cran.r-project.org/web/packages/policies.html
#
# see also: https://support.rstudio.com/hc/en-us/articles/200486518-Customizing-Package-Build-Options

if(.Platform$OS.type == "windows") {
  if(Sys.getenv("R_ARCH") == "/x64"){
    binR <- paste0("\"", Sys.getenv("R_HOME"), "/bin/x64/R.exe\"")
  } else{
    binR <- paste0("\"", Sys.getenv("R_HOME"), "/bin/R.exe\"")
  }
} else {
  binR <- "R"
}

version <- read.dcf("./DESCRIPTION")[, "Version"]

# build .tar.gz (with vignettes etc.) and check it afterwards
# the resulting postgwas.tar.gz can be sent to CRAN
setwd(paste0(getwd(), "/..")) # do not pollute our source directory
shell(paste0(binR, " CMD build postgwas"))  #takes very long time!
shell(paste0(binR, " CMD build postgwas --no-build-vignettes"))


# OPTIONAL: install and create *.zip windows binary package
shell(paste0(binR, " CMD INSTALL --build postgwas_", version, ".tar.gz"))

# run CRAN checks
shell(paste0(binR, " CMD check --as-cran postgwas_", version, ".tar.gz"))
setwd(paste0(getwd(), "/postgwas"))


# check reverse dependencies to inform maintainers of packages with reverse dependencies
source("http://developer.r-project.org/CRAN/Scripts/depends.R")
reverse_dependencies_with_maintainers("postgwas")

# release package to CRAN with the help of devtools function
devtools::release()

# check results of CRAN processing 48h later under:
# http://cran.r-project.org/web/checks/check_results_postgwas.html
