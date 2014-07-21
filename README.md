
# postgwas #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version 1.11-2
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

## Installation ##

### Install from CRAN ###

```R
install.packages("postgwas")
library(postgwas)
```
See [CRAN package repository](http://cran.r-project.org/web/packages/postgwas/index.html)


### Install from repository ### 

1. Install `devtools` via CRAN:
```R
install.packages(c("devtools", "rstudioapi"))
```
You might also want to [update devtools](https://github.com/hadley/devtools#updating-to-the-latest-version-of-devtools) to the newest developer version available.

2. Install `postgwas` dependencies not on CRAN from Bioconductor:

```R
biocPkgs <- c("biomaRt", "grid", "AnnotationDbi", "org.Hs.eg.db", "Rgraphviz", "RBGL")
source("http://bioconductor.org/biocLite.R")
biocLite(biocPkgs)
```

3. Install (if not already done) a working development environment:
    * **Linux**: Install a compiler for your distribution. For instance for Ubuntu this would be `sudo apt-get install r-base-dev`. Further instructions can be found at [CRAN](http://cran.r-project.org/bin/linux).
    * **Windows**: You will need the [Rtools](http://cran.r-project.org/bin/windows/Rtools/) installed.
    * **Mac**: You can find the latest Xcode which comes with C compilers in the Mac App Store.

4. Install `postgwas` from Bitbucket via `devtools`:

```R
devtools::install_bitbucket("postgwas")
```

## Contribution ##

You are free to contribute!

### Who do I talk to? ###

* Repo owner 
* Other community or team contact
