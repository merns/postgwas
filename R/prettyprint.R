prettyGene <- function(gene) {
  gene <- gsub("\\s", "", gene)
  if( is.null(gene["genename"]) || is.na(gene["genename"]) )
    return(NA)
  else
    paste(gene["genename"], "[", paste(gene["start"], gene["end"], sep = ":"), "]")
}

pasteNAasEmpty <- function(x, y, sep) {
  if(length(x) > 1 | length(y) > 1) {
    mapply(pasteNAasEmpty, x, y, sep)
  } else {
    nax <- is.na(x) | x == "NA"
    nay <- is.na(y) | y == "NA"
    if(nax & nay)
      return(NA)
    else if(nax)
      return(y)
    else if(nay)
      return(x)    
    else
      return(paste(x, y, sep = sep))
  }
}

# cluster is a list of data frames of gene entries
prettyCluster <- function(cluster) {
  if(length(cluster) < 1)
    return(NA)
  prettyHead   <- apply(cluster[[1]], 1, prettyGene)
  if(length(cluster) == 1)
    return(prettyHead)
  else
    pasteNAasEmpty(prettyHead, prettyCluster(cluster[-1]), " | ")
}


