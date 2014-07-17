list2df <- function(lst) {
	# Concatenates list elements, where elements have to be data frames that match in columns. Their rows will be concatenated.
	# when list elements are vectors, a matrix will be returned.
	# 
	# Args:
	#  lst: A list with nesting level of one, containing data frames or vectors.
	#
	# Value:
	#  a data frame or matrix, with input list elements concatenated row-wise.
	
	med <- floor(length(lst)/2)
	if(length(lst) > 1)
		rbind(list2df(lst[1:med]), list2df(lst[(med+1):length(lst)]))
	else
		return(rbind(lst[[1]])) # the rbind will convert vectors to data frames, but leaves data frames unchanged
}


na.set <- function(log.vec, set = FALSE) {
  # in a vector of booleans log.vec, set all NAs to TRUE or FALSE
  
  log.vec[is.na(log.vec)] <- set
  return(log.vec)
}

vectorElements <- function(df) {
  # cast all columns of a data frame to vector class
  
  if(!is.data.frame(df)) stop("Parameter is not a data frame")
  for(idx in 1:ncol(df)) df[, idx] <- as.vector(df[, idx])
  return(df)
}

lexsortColumns <- function(df, col1, col2) {
  # changes entries between col1 and col2 so that col1 < col2
  to.col1           <- df[, col1] > df[, col2]
  swap              <- df[to.col1, col1]
  df[to.col1, col1] <- df[to.col1, col2]
  df[to.col1, col2] <- swap
  return(df)
}

pruneVec <- function(x, max) {
	# Removes (evenly distributed) elements from vector x so that 'max' elements remain
	# 
	# Args:
	#  x:   A vector which is cut to 'max' number of elements
	#  max: the number of elements in the returned vector
	#
	# Value:
	#  A vector
	
	if(length(x) > max) {
		idx.length <- length(x)
		idx.false <- length(x) - max
		prune.idx <- rep(TRUE, idx.length)
		prune.idx[as.integer(1:idx.false * (idx.length / idx.false))] <- FALSE
		return(x[prune.idx])
	} else {
		return(x)
	}
}

makeBlocks <- function(x, block.count = NULL, block.size = sqrt(length(x))) {
	# Divides the data (list or vector) into blocks.
	# Division can be done controlling the number of resulting blocks or the (maximum) size of blocks. 
	#
	# Args:
	#  x:           A vector or list of data
	#  block.count: The number of blocks to divide the data into. Overrides block.size modus when not NULL.
	#  block.size:  Will divide the data into the smalles number of blocks so that each block has at most block.size elements. Is not applied when block.count is set.
	#
	# Value:
	#  A list, where each element contains a reduced number of elements from the input vector.
	#  Applying unlist to the returned list will give all elements (unordered) from the input data
		
	if(is.null(block.count) && is.null(block.size))
		stop("Either block.count or block.size has to be specified!")
	
	if(!is.null(block.count)) {
    idx <- rep(1:block.count, length.out = length(x))
    return(split(x, idx))
  } else {
    idx <- ceiling(seq_along(x) / block.size)
    return(split(x, idx))
  }
	
}


# takes a character vector 'sentences'
# when the number of characters at each element is longer than maxchar, 
# the spaces in the sentence will be counted and a newline inserted before the median space position
# returns the converted sentences vector
splitSentence <- function(sentences, maxchar = 40) {
  if(length(sentences) < 1 || !is.character(sentences))
    return(sentences)
  # linebreak long term labels
  needs.split <- nchar(sentences) > maxchar
  splithere <- sapply(
      gregexpr(" |-", sentences), # gives a list of vectors of space positions in each term name 
      function(spacepos) {
        # return the median space position, we will insert a newline there
        if(length(spacepos) > 0)
          spacepos[ceiling((length(spacepos)+1)/2)]
        else
          NA
      }
  )
  needs.split[is.na(splithere)] <- FALSE
  sentences[needs.split] <- paste(
      substr(sentences[needs.split], 1, splithere[needs.split]),
      substr(sentences[needs.split], splithere[needs.split] +1, nchar(sentences[needs.split])),
      sep = "\n"
  )
  return(sentences)
} 


nextFilename <- function(basename, extension, startcount = 1) {
	# Looks for an unoccupied filename in the current directory, 
	# starting with the name [basename][startcount].[extension] 
	# then increasing startcount until a free filename is found.
		
	if(!is.numeric(startcount))
		stop("Filename startcount has to be numeric")
	
	if(file.exists(paste(basename, startcount, ".", extension, sep = "")))
		return(nextFilename(basename, extension, startcount = startcount+1))
	else
		return(paste(basename, startcount, ".", extension, sep = ""))
	
}


getDensity <- function(elements, win, pos.colname = "mean") {
  # uses a sliding window approach to check how many elements (max) can be found within such a window
  elements <- elements[order(elements[, pos.colname]), ]
  max(sapply( 
      1:nrow(elements),      # for each element, return the number of elements that are within 'win' range
      function(idx.row) {
        sum(
          elements[, pos.colname] >= elements[idx.row, pos.colname] & 
          elements[, pos.colname] < elements[idx.row, pos.colname] + win
        )
      }
  ))
}

