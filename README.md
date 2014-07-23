# postgwas #

Facilitates annotation of genes to SNPs using proximity or LD
information, creates regional and manhattan plots and contains an
interaction network analysis tool for GWAS result data. Special features
cover subphenotype (intermediate phenotype) comparison and rare variant
display.

* See the [Vignette](src/inst/doc/postgwas.pdf) for further information and a lot of examples. 
* There is also a [publication](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0071775) covering `postgwas`.

* * *

## Installation ##

### Option 1: Install from CRAN (1 minute) ###

    install.packages("postgwas")
    library(postgwas)

See also the [CRAN package repository](http://cran.r-project.org/web/packages/postgwas/index.html).


### Option 2: Install from Bitbucket (5 minutes) ###

1. Install (if you haven't already) a working development environment:
    * **Linux**: Install a compiler for your distribution. For instance for Ubuntu this would be `sudo apt-get install r-base-dev`. Further instructions can be found at [CRAN](http://cran.r-project.org/bin/linux).
    * **Windows**: You will need to have the [Rtools](http://cran.r-project.org/bin/windows/Rtools/) installed.
    * **Mac**: You can find the latest Xcode.app containing a C compiler in the Mac App Store.

2. Install `devtools` via CRAN:

        install.packages(c("devtools", "rstudioapi"))
        
    (Note: You might also want to [update devtools](https://github.com/hadley/devtools#updating-to-the-latest-version-of-devtools) to the newest version available.)
    
3. Install dependencies of `postgwas` that are not on CRAN from Bioconductor:

        biocPkgs <- c("biomaRt", "AnnotationDbi", "org.Hs.eg.db", "Rgraphviz", "RBGL")
        source("http://bioconductor.org/biocLite.R")
        biocLite(biocPkgs)

4. Install `postgwas` from Bitbucket via `devtools`:

        devtools::install_bitbucket("postgwas")


* * *

## Contribution ##

You are welcome to contribute!

Just contact one of the Repo owners [Marko Ernsting](https://bitbucket.org/merns), [Milan Hiersche](mailto:mhiersche@gmx.de) or [Frank Rühle](https://www.researchgate.net/profile/Frank_Ruehle).

* * *

## Citation ##

If you use the package for actual research, please cite the following [PlosOne publication](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0071775):

> Hiersche, M., Ruehle, F., & Stoll, M. (2013). Postgwas: Advanced GWAS 
> Interpretation in R. PloS one, 8(8), e71775. doi:10.1371/journal.pone.0071775


Use the following BibTex entry or download citation information from [here](http://www.plosone.org/article/citationList.action?articleURI=info%3Adoi%2F10.1371%2Fjournal.pone.0071775).
```Latex
@article{10.1371/journal.pone.0071775,
    author = {Hiersche, , Milan AND Rühle, , Frank AND Stoll, , Monika},
    journal = {PLoS ONE},
    publisher = {Public Library of Science},
    title = {Postgwas: Advanced GWAS Interpretation in R},
    year = {2013},
    month = {08},
    volume = {8},
    url = {http://dx.doi.org/10.1371%2Fjournal.pone.0071775},
    pages = {e71775},
    abstract = {We present a comprehensive toolkit for post-processing, visualization and advanced analysis of GWAS results. In the spirit of comparable tools for gene-expression analysis, we attempt to unify and simplify several procedures that are essential for the interpretation of GWAS results. This includes the generation of advanced Manhattan and regional association plots including rare variant display as well as novel interaction network analysis tools for the investigation of systems-biology aspects. Our package supports virtually all model organisms and represents the first cohesive implementation of such tools for the popular language R. Previous software of that range is dispersed over a wide range of platforms and mostly not adaptable for custom work pipelines. We demonstrate the utility of this package by providing an example workflow on a publicly available dataset.},
    number = {8},
    doi = {10.1371/journal.pone.0071775}
}        
```