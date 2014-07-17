# initialize empty buffer variable list
# the variables in the environment postgwasBuffer always have to exist (and do so by implementation)
# they can be NULL which means 'not set'
postgwasBuffer <- new.env()
GOSimEnv <- new.env()
assign("snps", NULL, envir = postgwasBuffer)
assign("genes", NULL, envir = postgwasBuffer)
assign("genes.regionalplot", NULL, envir = postgwasBuffer)
assign("exons.regionalplot", NULL, envir = postgwasBuffer)
assign("ld.regionalplot", NULL, envir = postgwasBuffer)
assign("goterms", NULL, envir = postgwasBuffer)
