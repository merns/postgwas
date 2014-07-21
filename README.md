# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
  [CRAN package repositors](http://cran.r-project.org/web/packages/postgwas/index.html)
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact

```
#!S
#
devtools::install_bitbucket("postgwas", username="merns", password="")
version <- read.dcf("./DESCRIPTION")[, "Version"]
```

```
#!R
#
devtools::install_bitbucket("postgwas", username="merns", password="")
version <- read.dcf("./DESCRIPTION")[, "Version"]
```


~~~~
#
devtools::install_bitbucket("postgwas", username="merns", password="")
version <- read.dcf("./DESCRIPTION")[, "Version"]
~~~~
