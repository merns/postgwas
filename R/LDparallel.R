LDparallel <- function(g, num.processes, ...) {
	
	require(parallel)
  if(!'genetics' %in% installed.packages()[, 'Package'])
    stop("getGenotypes: Package 'genetics' has to be installed for parallelized LD calculation")
  else 
    require(genetics)
  
	if(missing(num.processes))
		num.processes <- detectCores()
	
	# check parameter validity
	gvars <- sapply( g, function(x) (is.genotype(x) && nallele(x)==2) )
	if(any(gvars==FALSE)){
		warning("Non-genotype variables or genotype variables ",
				"with more or less than two alleles detected. ",
				"These variables will be omitted: ",                
				paste( colnames(g)[!gvars] , collapse=", " )
		)
		g <- g[,gvars]
	}
	
	# matrix with to-be calculated LD values, dimnames == input data dimnames
	ld.res <- matrix(nrow = length(g), ncol = length(g))
	rownames(ld.res) <- colnames(g)
	colnames(ld.res) <- colnames(g)
	
	# Parallelisiere Berechnung fuer jede Zeile der Matrix, da sonst ncol(g)*ncol(g) Prozesse gestartet werden
	# Verteile arbeit auf CPUs, Jeder Kern erhaelt (ncol(g)-1) / cores Zeilen, wobei die Aufteilung der Zeilen in cores spruengen passiert, 
	# da die zeilenberechnung immer weniger zeilen enthaelt
	# Dies ist die sinnvollste variante, da nur die Haelfte der Matrix berechnet werden muss und die Anzahl der zu berechnenden Werte mit steigender Zeilennummer jeweils sinkt
	
	# generate all pairs of snps (as index pairs referencing ld.res as well as g columns)
	idx.pairs <- merge(1:length(g), 1:length(g))
	# remove symmetric indices (e.g. keep only one pair of (2,4) and (4,2))
	idx.pairs <- idx.pairs[idx.pairs$x > idx.pairs$y, ]
	idx.pairs.list <- apply(idx.pairs, 1, list)
	
	res <- mclapply(
			idx.pairs.list, 
			function(idx.pair) {
				x <- idx.pair[[1]]["x"]
				y <- idx.pair[[1]]["y"]
				return(list(x = x, y = y, ld = LDparallel.genotype(g[[x]], g[[y]])$"R^2"))
			}, 
			mc.cores = num.processes
	)
	
	for(pair in res)
		ld.res[pair$y, pair$x] <- pair$ld
	
	
	return(ld.res)
	
}


LDparallel.genotype <- function(g1,g2,...) {
	
	if(is.haplotype(g1) || is.haplotype(g2))
		stop("Haplotype options are not yet supported.")
	
	if(nallele(g1)!=2 || nallele(g2)!=2)
		stop("This function currently only supports 2-allele genotypes.")
	
	prop.A <- summary(g1)$allele.freq[,2]
	prop.B <- summary(g2)$allele.freq[,2]
	
	major.A <- names(prop.A)[which.max(prop.A)]
	major.B <- names(prop.B)[which.max(prop.B)]
	pA <- max(prop.A, na.rm=TRUE)
	pB <- max(prop.B, na.rm=TRUE)
	pa <- 1-pA
	pb <- 1-pB
	
	Dmin <- max(-pA*pB, -pa*pb)
	pmin <- pA*pB + Dmin;
	
	Dmax <- min(pA*pb, pB*pa);
	pmax <- pA*pB + Dmax;
	
	counts <- table(
			allele.count(g1, major.A),
			allele.count(g2, major.B)
	)
	
	n3x3 <- matrix(0, nrow=3, ncol=3)
	colnames(n3x3) <- rownames(n3x3) <- 0:2
	
	# ensure the matrix is 3x3, with highest frequency values in upper left
	for(i in rownames(counts))
		for(j in colnames(counts))
			n3x3[3-as.numeric(i),3-as.numeric(j)] <- counts[i,j]
	
	
	loglik <- function(pAB,...)
	{
		(2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
				(2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(pA-pAB) +
				(2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(pB-pAB) +
				(2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-pA-pB+pAB) + 
				n3x3[2,2]*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB))
	}
	
	# SAS code uses:
	#
	# s <- seq(pmin+0.0001,pmax-0.0001,by=0.0001)
	# lldmx <- loglik(s)
	# maxi <- which.max(lldmx)
	# pAB <- s[maxi]
	
	# but this should be faster:
	solution <- optimize(
			loglik,
			lower=pmin+.Machine$double.eps,
			upper=pmax-.Machine$double.eps,
			maximum=TRUE
	)
	pAB <- solution$maximum
	
	estD <- pAB - pA*pB
	if (estD>0)  
		estDp <- estD / Dmax
	else
		estDp <- estD / Dmin
	
	n <-  sum(n3x3)
	
	corr <- estD / sqrt( pA * pB * pa * pb )
	
	# dchi <- (2*n*estD^2)/(pA * pa * pB* pb)
	# dpval <- 1 - pchisq(dchi,1)
	
	retval <- list(
			call=match.call(),
			"D"=estD,
			"D'"=estDp,
			"r" = corr,
			"R^2" = corr^2,
			"n"=n
	    # "X^2"=dchi,
	    # "P-value"=dpval
	)
	
	class(retval) <- "LD"
	retval
	
}
