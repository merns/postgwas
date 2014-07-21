
# update devtools ==========================================================

install.packages(c("devtools", "rstudioapi"))
devtools::build_github_devtools()
#### Restart R before continuing ####
install.packages("./devtools.zip", repos = NULL)
unlink("./devtools.zip")


# Development with "devtools" ===============================================

devtools::dev_mode(TRUE)

# loading the package
devtools::load_all()
devtools::test()


devtools::dev_mode(FALSE)






