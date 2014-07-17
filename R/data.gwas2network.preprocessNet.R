
# assumes a variable 'network' in environment postgwasBuffer (a data frame) 

# value: n/a, transforms 'network' as a side effect in postgwasBuffer envir:
# it then has uppercase geneids/names in col 1 and 2 with the lexicographic smaller names in col 1, 
# a weight column
# and all columns are character vectors (not factors).
data.gwas2network.preprocessNet <- function() {

  message("Preprocessing network data... ", appendLF = FALSE)
  net <- get("network", envir = postgwasBuffer)
  
  # make sure weight column exists and does not contain NA (set to 1 instead)
  if(is.null(net$weight))
    net$weight <- 1
  net$weight[is.na(net$weight)] <- 1
  net <- na.omit(net[, c(1, 2, which(colnames(net) %in% c("weight", "label")))])
  if(nrow(net) < 1)
    stop("Network has no edges!\n")
  
  net[, 1] <- toupper(as.character(as.vector(net[, 1])))
  net[, 2] <- toupper(as.character(as.vector(net[, 2])))
  if(ncol(net) > 2)
    net[, 3] <- as.character(as.vector(net[, 3]))
  if(ncol(net) > 3)
    net[, 4] <- as.character(as.vector(net[, 4]))
  
  # put lexicographically smaller gene names in column 1 (all gene names in col2 will be larger than at the corresponding row in col1)
  from2to1 <- net[, 1] > net[, 2]
  swap <- net[from2to1, 1]
  net[from2to1, 1] <- net[from2to1, 2]
  net[from2to1, 2] <- swap
  net[, 1] <- as.character(as.vector(net[, 1]))
  net[, 2] <- as.character(as.vector(net[, 2]))
  
  assign("network", net, envir = postgwasBuffer)
  message("done.")
}