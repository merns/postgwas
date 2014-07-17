## The following code has been extracted from the GOSim package by Holger Fröhlich et. al. and adapted for parallelization.
## Holger Froehlich, N. Speer, A. Poustka, Tim Beissbarth (2007). GOSim - An R-Package for Computation of Information Theoretic GO Similarities Between Terms and Gene Products. BMC Bioinformatics, 8:166

# req in setOntology
getAncestors<-function(){
  ontology<-get("ontology",envir=GOSimEnv)
  if(ontology == "BP")
    res<-AnnotationDbi::as.list(GOBPANCESTOR)
  else if(ontology == "MF")
    res<-AnnotationDbi::as.list(GOMFANCESTOR)
  else if(ontology == "CC")
    res<-AnnotationDbi::as.list(GOCCANCESTOR)
  else
    stop(paste("ontology", ontology, "not known!"))
  return(res)
}

# req in setOntology
getParents<-function(){
  ontology<-get("ontology",envir=GOSimEnv)
  if(ontology == "BP")
    res<-AnnotationDbi::as.list(GOBPPARENTS)
  else if(ontology == "MF")
    res<-AnnotationDbi::as.list(GOMFPARENTS)
  else if(ontology == "CC")
    res<-AnnotationDbi::as.list(GOCCPARENTS)
  else
    stop(paste("ontology", ontology, "not known!"))
  return(res)
}

# req in setOntology
getChildren<-function(){
  ontology<-get("ontology",envir=GOSimEnv)
  if(ontology == "BP")
    res<-AnnotationDbi::as.list(GOBPCHILDREN)
  else if(ontology == "MF")
    res<-AnnotationDbi::as.list(GOMFCHILDREN)
  else if(ontology == "CC")
    res<-AnnotationDbi::as.list(GOCCCHILDREN)
  else
    stop(paste("ontology", ontology, "not known!"))
  return(res)
}

# req in calcICs
getOffsprings<-function(){
  ontology<-get("ontology",envir=GOSimEnv)
  if(ontology == "BP")
    res<-AnnotationDbi::as.list(GOBPOFFSPRING)
  else if(ontology == "MF")
    res<-AnnotationDbi::as.list(GOMFOFFSPRING)
  else if(ontology == "CC")
    res<-AnnotationDbi::as.list(GOCCOFFSPRING)
  else
    stop(paste("ontology", ontology, "not known!"))
  return(res)
}

# req in setOntology
calcICs<-function(cores = 1){	
#  if(!require(annotate))
#    stop("Package annotate is required for function calcICs")
  evidences<-get("evidences", envir=GOSimEnv)
  ontology<-get("ontology",envir=GOSimEnv)
  organism = get("organism", envir=GOSimEnv)
  message(paste("calculating information contents for ontology ", ontology, " using evidence codes '", paste(evidences,collapse=", "), "' (", organism,") ...",sep=""))	
  ids<- AnnotationDbi::toTable(GOTERM)	 
  ids = unique(ids[ids[,"Ontology"] == ontology,"go_id"]) # these are all GO terms, which belong to the corrseponding category
  offspring <- getOffsprings()	
  gomap <- get("gomap",envir=GOSimEnv)		
  goterms<-unlist(sapply(gomap, function(x) names(x)), use.names=FALSE) # all GO terms appearing in an annotation
  goterms<-goterms[goterms %in% ids] # this is to ensure that we only get GO terms mapping to the given ontology
  tab<-table(goterms)
  na<-names(tab)			
  s<-setdiff(ids, na)  #ensure that GO terms NOT appearing in the annotation have 0 frequency
  m<-double(length(s))
  names(m)<-s
  tab<-as.vector(tab)
  names(tab)<-na
  tab<-c(tab,m)
  if(cores <= 1) {
    ta<-sapply(ids,function(x){ t=tab[unlist(offspring[x])]; tab[x]+sum(t[!is.na(t)])})
  } else {
    ta<-simplify2array(mclapply(ids,function(x){ t=tab[unlist(offspring[x])]; tab[x]+sum(t[!is.na(t)])}, mc.cores = cores))    
  }
  names(ta)<-ids	
  IC <- -log(ta/sum(tab))	# ACHTUNG: hier muß tab und nicht ta stehen: Die Vorkommenshäufigkeit des Wurzelknotens ist die Summe der Vorkommenshäufigkeiten aller Knoten OHNE Aufsummieren der Kinder!
# # 	IC[IC == Inf] = 0 # WRONG: GO terms which are not annotated have Inf information content (NOT 0: They cannot be treated like root!!!)
  message("done")
  return(IC)
}

setEvidenceLevel<-function(evidences="all", organism=org.Hs.egORGANISM, gomap=org.Hs.egGO, cores = 1){			
  message(paste("-> retrieving GO information for all available genes for organism '", organism, "' in GO database", sep=""))
  assign("evidences", evidences, envir=GOSimEnv)	
# 	gomap<-as.list(GOENTREZID2GO)		
  if(is(gomap, "Bimap")){		
    mapped_genes <- mappedkeys(gomap)	
    gomap = AnnotationDbi::as.list(gomap[mapped_genes])
  }
  else if(!is(gomap, "list"))
    stop("gomap argument should be a nested list (see manual pages)!")
  message(paste("-> filtering GO terms according to evidence levels '", evidences, "'",sep=""))
  if((length(evidences) > 1) || (evidences!="all")){		
    gomap<-sapply(gomap,function(x) sapply(x,function(y) c(y$Evidence %in% evidences, y$Ontology)))
    gomap<-sapply(gomap, function(x) x[2,x[1,]=="TRUE"])
    gomap<-gomap[sapply(gomap,length) >0]
  }
  assign("gomap", gomap, envir=GOSimEnv)
  assign("organism", organism, envir=GOSimEnv)	
}

setOntology<-function(ont="BP", loadIC=FALSE, cores = 1){
  assign("ontology", ont, envir=GOSimEnv)
  if(loadIC){
    if(!exists("IC",envir=GOSimEnv)) 
      stop('Argument loadIC is TRUE, but IC object does not exist in GOSimEnv')
    IC<-get("IC",envir=GOSimEnv)
    IC["all"]=0
  } else {
    IC <- calcICs(cores = cores)
    assign("IC", IC, envir=GOSimEnv)
  }
  assign("ancestor", getAncestors(), envir=GOSimEnv) 	
  assign("children", getChildren(), envir=GOSimEnv)
  assign("parents", getParents(), envir=GOSimEnv) 
  children<-get("children",envir=GOSimEnv) 	
  parents<-get("parents",envir=GOSimEnv) 	
  assign("nchildren", sapply(children,length) , envir=GOSimEnv)
  assign("nparents", sapply(parents,length), envir=GOSimEnv)
  nchildren<-get("nchildren",envir=GOSimEnv) 	
  nparents<-get("nparents",envir=GOSimEnv) 	
  assign("Eavg", (sum(nchildren) + sum(nparents))/ length(union(names(nchildren), names(nparents))), envir<-GOSimEnv)
  assign("alphaParam", 0.5, envir=GOSimEnv)
  assign("betaParam", 0.5, envir=GOSimEnv)  	    	
}

initialize <- function(evidences = "all", organism = "human", ontology = "BP", gomap=org.Hs.egGO, cores = 1) {
  if(!exists('GOSimEnv'))
    stop('There is no environment "GOSimEnv" defined, cannot initialize modified GOSim code')
  GOSimEnv$organism <- NULL
  GOSimEnv$ontology <- NULL
  require(GO.db)
}







# req by get gene sim
# filter out genes not mapping to the category in question	
filterGO<-function(genelist){	
  cat("filtering out genes not mapping to the currently set GO category ...")
  if(!exists("GOSimEnv")) initialize()	
  IC<-get("IC", envir=GOSimEnv)
  ids<-names(IC[IC != Inf]) # only consider GO terms with some known annotation   
  gomap<-get("gomap",envir=GOSimEnv)
  k<-1
  allgenes<-list()
  for(i in 1:length(genelist)){		
    annoi<-gomap[[match(as.character(genelist[i]), names(gomap))]]
    annoi<-intersect(names(annoi), as.character(ids))				
    if(length(annoi) > 0){
      allgenes[[k]]<-list(annotation=annoi,genename=as.character(genelist[i]))
      k<-k + 1
    }
  }	
  cat(" ===> list of ", length(genelist), "genes reduced to ", length(allgenes), "\n")
  allgenes
}



# precompute term similarities for all pairs of GO terms belonging to the annotated genelists x (and y)
precomputeTermSims<-function(x, y=NULL, similarityTerm="JiangConrath", verbose=FALSE){ 	
  if(verbose)
    message("precomputing term similarities ...")	
  gotermsx<-as.vector(unique(unlist(sapply(x, function(xx) xx$annotation))))
  if(!is.null(y)){  		 
    gotermsy<-as.vector(unique(unlist(sapply(y, function(xx) xx$annotation))))
    if(similarityTerm %in% c("diffKernelgraphLapl", "diffKernelLLE", "diffKernelpower", "diffKernelexpm")){
      STerm = getTermSim(c(gotermsx, gotermsy), method=similarityTerm)
      STerm = STerm[gotermsx, gotermsy]
      return(STerm)
    }
    STerm<-matrix(0, nrow=length(gotermsx), ncol=length(gotermsy))
    rownames(STerm)=gotermsx
    colnames(STerm)=gotermsy
    for(i in 1:length(gotermsx)){						
      for(j in 1:length(gotermsy)){          
        STerm[i,j]<-calcTermSim(gotermsx[i],gotermsy[j], similarityTerm, verbose)
      }			
    }		
  } else{
    if(similarityTerm %in% c("diffKernelgraphLapl", "diffKernelLLE", "diffKernelpower", "diffKernelexpm")){
      STerm = getTermSim(gotermsx, method=similarityTerm)			
      return(STerm)
    }
    STerm<-matrix(0, nrow=length(gotermsx), ncol=length(gotermsx))
    rownames(STerm)<-gotermsx
    colnames(STerm)<-gotermsx
    for(i in 1:length(gotermsx)){
      STerm[i,i]<-calcTermSim(gotermsx[i], gotermsx[i], similarityTerm, verbose)
      if(i > 1){
        for(j in 1:(i-1)){					
          STerm[i,j]<-calcTermSim(gotermsx[i], gotermsx[j], similarityTerm, verbose)
          STerm[j,i]<-STerm[i,j]
        }
      }
    }
  }	
  STerm
}



# basic term similarity for index i and j in GOIDs ids
calcTermSim<-function(ids, i, j, method="JiangConrath"){	
  calcTermSim(ids[i],ids[j], method)
}

# basic term similarity between term1 and term2
calcTermSim<-function(term1, term2, method="JiangConrath", verbose=FALSE){
  if(!exists("GOSimEnv")) initialize()
  IC<-get("IC", envir=GOSimEnv)
  if(verbose)
    message(paste("Terms:",term1,",",term2,"( method:",method,")"))	
  if(method== "Resnik")
    return(IC[getMinimumSubsumer(term1,term2)] / max(IC[IC != Inf]))   
  else if(method == "JiangConrath")
    return(1 - min(1, -2*IC[getMinimumSubsumer(term1,term2)] + IC[term1] + IC[term2]) )	
  else if(method == "Lin"){
    res = 2*IC[getMinimumSubsumer(term1,term2)]/(IC[term1]+IC[term2])
    return(ifelse(is.na(res), 0, res))
  }
  else if(method== "CoutoEnriched")
    return(getEnrichedSim(term1, term2))
  else if(method == "CoutoResnik")  
    return(getDisjCommAncSim(term1,term2, "Resnik"))
  else if(method == "CoutoJiangConrath")  
    return(getDisjCommAncSim(term1,term2, "JiangConrath"))
  else if(method == "CoutoLin"){
    res = getDisjCommAncSim(term1,term2, "Lin")
    return(ifelse(is.na(res), 0, res))
  }
  else if(method == "simIC"){ # Li et al.
    MICA = getMinimumSubsumer(term1,term2)
    res = 2*IC[MICA]/(IC[term1]+IC[term2]) * (1 - 1/(1 + IC[MICA]))
    return(ifelse(is.na(res), 0, res))
  }
  else if(method == "GIC") # graph information content
    return(getGIC(term1, term2))
  else if(method == "relevance"){ # Schlicker et al.
    MICA = getMinimumSubsumer(term1,term2)
    res = (2*IC[MICA]/(IC[term1]+IC[term2]))*(1 - exp(-IC[MICA]))
    return(ifelse(is.na(res), 0, res))
  }
  else if(method == "diffKernel"){
    K = mget("K", envir=GOSimEnv, ifnotfound=list(function(x) stop(paste("Diffusion kernel not loaded!\nPlease invoke load.diffusion.kernel().", sep = ""), call. = FALSE)))$K
    return(K[term1, term2])
  }
  else
    stop(paste("calcTermSim: Unknown term similarity",method))
}


getMinimumSubsumer<-function(term1, term2){	
  if(!exists("GOSimEnv")) initialize()
  ancestor<-get("ancestor",envir=GOSimEnv)
  if(term1 == term2){
    ms<-term1	
    return(ms)
  }
  an1<-unlist(ancestor[names(ancestor) == term1])
  an2<-unlist(ancestor[names(ancestor) == term2])
  case1<-which(an2 == term1)  # term1 is the ms of term2
  case2<-which(an1 == term2) # term2 is the ms of term1	  
  if(length(case1) > 0){
    ms<-term1	
  } else if(length(case2) > 0) {
    ms<-term2	
  } else {
    # find common ancestor with maximal information content
    anall<-intersect(an1, an2) 
    IC<-get("IC", envir=GOSimEnv)		
    ms<-anall[which.max(IC[anall])]
  }	
  if(is.null(ms) | length(ms) == 0)
    ms <- "NA"	
  ms
}



# compute gene similarity for a pair of genes having GO terms anno1 and anno2
getGSim<-function(anno1, anno2, similarity="funSimMax", similarityTerm="JiangConrath", STerm=NULL, avg=FALSE, verbose=FALSE){	  
  if(length(anno1) <= length(anno2)){
    a1<-anno1
    a2<-anno2
    swap<-FALSE
  }
  else{
    a1<-anno2
    a2<-anno1
    swap<-TRUE
  }		
  if(!is.null(STerm)){ # use precomputed similarity values
    if(!swap)
      ker<-STerm[a1,a2]
    else
      ker<-STerm[a2,a1]		
    if(length(a1) == 1)
      ker<-t(as.matrix(ker))		
    if(is.null(ker) || is.null(nrow(ker))){
      warning(paste("No GO information for",a1,a2,". Similarity set to NaN."))		
      return(NaN)
    }	
    if(nrow(ker) > ncol(ker))
      ker<-t(ker)
  }
  else{ 
    if(similarity %in% c("dot"))
      return(getWeightedDotSim(a1, a2))			
    else{
      # calculate term similarity
      ker<-matrix(0,nrow=length(a1),ncol=length(a2))	
      for(i in 1:length(a1)){
        for(j in 1:length(a2))
          ker[i,j]<-calcTermSim(a1[i],a2[j], similarityTerm, verbose)		
      }
    }
  }  
  if(length(a1)*length(a2) > 0){
    if(similarity == "OA"){				
      res<-.C("OAWrapper", ker, nrow(ker), ncol(ker), as.integer(1), ret=double(1))$ret
      if(avg)
        res = res/length(a2)	
      return(res)
    }
    else if(similarity == "max"){				
      return(max(ker))
    }
    else if(similarity == "mean"){				
      return(mean(ker))
    }  
    else if(similarity == "funSimAvg"){
      rowMax = mean(apply(ker,1,max))
      colMax = mean(apply(ker,2,max))
      return(0.5*(rowMax + colMax))
    }
    else if(similarity == "funSimMax"){
      rowMax = mean(apply(ker,1,max))
      colMax = mean(apply(ker,2,max))
      return(max(rowMax, colMax))
    }
    else if(similarity == "hausdorff"){
      rowMax = min(apply(ker,1,max))
      colMax = min(apply(ker,2,max))
      return(min(rowMax, colMax))	
    }	
    else
      stop(paste("getGSim: Unknown gene similarity",similarity,"!"))
  }
  else{	
    warning(paste("No GO information for",a1,a2,". Similarity set to NaN."))		
    return(NaN)
  }
}

# req for getGeneSim
normalize.kernel = function(Ker, kerself1=NULL, kerself2=NULL, method="none"){
  if(method != "none"){
    if(is.null(kerself1) || is.null(kerself2)){
      if(method == "sqrt"){ # result between -1 and 1
        Kd<-sqrt(diag(Ker) + 1e-10)
        Ker<-Ker/(Kd%*%t(Kd))			
      }
      else if(method == "Lin"){ # result: diagonal = 1		
        Kd = diag(Ker)
        Ker = 2*Ker / outer(Kd, Kd, "+")
      }
      else if(method == "Tanimoto"){ 
        Kd = diag(Ker)
        Ker = Ker / (outer(Kd, Kd, "+") - Ker)
      }			
#			else if(method == "variance")
#				Ker = Ker /(mean(diag(Ker)) - mean(Ker)) 
      else
        stop(paste("Unknown normalization method", method))
      diag(Ker) = 1
    }
    else{
      if(method == "sqrt")
        return(Ker / sqrt(kerself1 * kerself2))
      else if(method == "Lin")
        return(2*Ker / (kerself1 + kerself2))
      else if (method == "Tanimoto")
        return(Ker / (kerself1 + kerself2 - Ker))			
      else
        stop(paste("Unknown normalization method", method))
    }
  }
  Ker
}


getGeneSim<-function(genelist1, genelist2=NULL, similarity="funSimMax", similarityTerm="relevance", normalization=TRUE, method="sqrt", avg=(similarity=="OA"), verbose=FALSE, cores = 1){	
  genelist1 <- unique(genelist1)
  genelist2 <- unique(genelist2)
  if(length(genelist1) < 2 && is.null(genelist2))
    stop("Gene list should contain more than 2 elements!")
  allgenes<-filterGO(genelist1)	
  if(length(allgenes) > 1 && is.null(genelist2)){			
    if(!(similarity %in% c("dot")))
      STerm<-precomputeTermSims(x=allgenes, similarityTerm=similarityTerm, verbose=verbose) # precompute term similarities => speed up!		
    else
      STerm = NULL
    if(verbose)
      message(paste("Calculating similarity matrix with similarity measure",similarity))
    # computes the columns of the triangular part of matrix Ker and returns each row as a vector
    # this can be done in parallel
    tri.cols <- mclapply(
        seq_along(allgenes), 
        function(i) {
          tri.col <- numeric(length(allgenes)) # initialize new empty col 
          annoi <- allgenes[[i]]$annotation	
          for(j in 1:i) {
            annoj <- allgenes[[j]]$annotation
            tri.col[j] <- getGSim(annoi,annoj, similarity, similarityTerm, STerm=STerm, avg=avg, verbose)
          }
          return(tri.col)
        }, 
        mc.cores = cores
    )
    # populate Ker by setting rows and columns = rows (triangular matrix)
    Ker <- matrix(unlist(tri.cols), nrow=length(allgenes), ncol=length(allgenes))
    Ker[lower.tri(Ker)] <- t(Ker)[lower.tri(Ker)] # copies the upper triangle to the lower
    colnames(Ker) <- sapply(allgenes, function(x) x$genename)
    rownames(Ker) <- colnames(Ker)
    if(normalization){			
      Ker = normalize.kernel(Ker, method=method)
      if(any(Ker > 1, na.rm=T)) # this has been updated
        warning("Similarity matrix contains values > 1! This may happen with simlarity='funSimMax', if one gene's GO annotation is a complete subset of another gene's GO annotation.")
      Ker[Ker>1] = 1 # can happen with similarity funSimMax in cases where one GO annotation is subset of another one
    }			
  }
  else if(length(allgenes) > 0 && length(genelist2) > 0){
    allgenes2<-setdiff(filterGO(genelist2), allgenes)
    if(length(allgenes) < 1)
      stop('The supplied gene sets seem to be indentical')
    if(!(similarity %in% c("dot")))
      STerm<-precomputeTermSims(x=allgenes, y=allgenes2, similarityTerm=similarityTerm, verbose=verbose) # precompute term similarities => speed up!		
    else
      STerm = NULL
    if(verbose)
      message(paste("Calculating similarity matrix with similarity measure",similarity))
    Ker.list <- mclapply(
        seq_along(allgenes2), 
        function(i) {
          tri.col <- numeric(length(allgenes)) # initialize new empty col 
          annoi <- allgenes2[[i]]$annotation
          if(normalization)
            kerselfi = getGSim(annoi, annoi, similarity, similarityTerm, STerm=NULL, avg=avg, verbose)
          for(j in seq_along(allgenes)) {
            annoj <- allgenes[[j]]$annotation
            tri.col[j] <- getGSim(annoi,annoj, similarity, similarityTerm, STerm=STerm, avg=avg, verbose)
            if(normalization){
              kerselfj = getGSim(annoj, annoj, similarity, similarityTerm, STerm=NULL, avg=avg, verbose)
              tri.col[j] = normalize.kernel(tri.col[j], kerselfi, kerselfj, method=method)
            }
          }
          return(tri.col)
        }, 
        mc.cores = cores
    )
    Ker <- matrix(unlist(Ker.list), nrow=length(allgenes), ncol=length(allgenes2))
    colnames(Ker)<-sapply(allgenes2,function(x) x$genename)
    rownames(Ker)<-sapply(allgenes,function(x) x$genename)
    if(any(Ker > 1, na.rm=T))
      warning("Similarity matrix contains values > 1! This may happen with simlarity='funSimMax', if one gene's GO annotation is a complete subset of another gene's GO annotation.")
    Ker[Ker>1] = 1 # can happen with similarity funSimMax in cases where one GO annotation is subset of another one
  }
  else{
    if(length(allgenes) == 0)
      stop("No gene has GO information!")					
    else if(length(allgenes) == 1)
      stop(paste("Only gene",allgenes," has GO information!"))					
  }
  Ker
}






# compute FuSSiMeg enriched term similarity
getEnrichedSim<-function(term1, term2){   
  #if(!require(RBGL))
  #	stop("Package RBGL is required for function getDepthFactor")
  if(!exists("GOSimEnv")) initialize() 
  ms<-getMinimumSubsumer(term1,term2)
  IC<-get("IC", envir=GOSimEnv)	
  if(term1 != term2){    	    	
    G<-getGOGraph(c(term1,term2))
    if(term1 != ms){
      path1=sp.between(G,term1,ms)[[1]]$path # path to ms                
      len<-length(path1)
      delta1<-sum(sapply(path1[2:len],getDepthFactor,G)*sapply(path1[2:len],getDensityFactor)*(-diff(IC[path1])))
    }
    else
      delta1<-0
    if(term2 != ms){
      path2<-sp.between(G,term2,ms)[[1]]$path # path to ms    	
      len<-length(path2)
      delta2<-sum(sapply(path2[2:len],getDepthFactor,G)*sapply(path2[2:len],getDensityFactor)*(-diff(IC[path2])))
    }
    else
      delta2<-0
    delta<-delta1 + delta2		
    sim<-1 - min(1, delta)
  }
  else
    sim<-1 
  sim<-sim * IC[term1] * IC[term2]  # correction given in equation (11) of the FuSSiMeg paper
  sim[is.na(sim)] = 0
  names(sim)<-c()   
  sim
}


# get GraSM similarity of common disjunctive ancestors of two terms
getDisjCommAncSim<-function(term1, term2, method="JiangConrath"){
  if(!exists("GOSimEnv")) initialize()
  IC<-get("IC", envir=GOSimEnv)
  djca<-getDisjCommAnc(term1, term2)	
  ICdjca<-IC[djca]	
  ICdjca<-ICdjca[!is.na(ICdjca)]							
  ICshare<-mean(ICdjca)		
  if(method == "JiangConrath")
    return(1-min(1,-2*ICshare + IC[term1] + IC[term2]))
  else if(method == "Resnik")
    return(ICshare)
  else if(method == "Lin")
    return(2*ICshare/(IC[term1]+IC[term2]))
  else
    stop(paste("getDisjCommAnc: Unknown term similarity",method))
}

# graph information content similarity related to Tanimoto-Jacard index
getGIC = function(term1, term2){
  if(!exists("GOSimEnv")) initialize()	
  if(term1 == term2){
    return(1)
  }
  IC<-get("IC", envir=GOSimEnv)
  ancestor<-get("ancestor",envir=GOSimEnv)
  an1<-unlist(ancestor[names(ancestor) == term1])
  an2<-unlist(ancestor[names(ancestor) == term2])
  ancommon = intersect(an1, an2)
  anunion = union(an1, an2)
  return(sum(IC[ancommon]) / sum(IC[anunion]))
}


# graph information content similarity related to Tanimoto-Jacard index
getGIC = function(term1, term2){
  if(!exists("GOSimEnv")) initialize()	
  if(term1 == term2){
    return(1)
  }
  IC<-get("IC", envir=GOSimEnv)
  ancestor<-get("ancestor",envir=GOSimEnv)
  an1<-unlist(ancestor[names(ancestor) == term1])
  an2<-unlist(ancestor[names(ancestor) == term2])
  ancommon = intersect(an1, an2)
  anunion = union(an1, an2)
  return(sum(IC[ancommon]) / sum(IC[anunion]))
}

getWeightedDotSim <- function(anno1, anno2){
  v1 = getGeneFeatures.internal(anno1)
  v2 = getGeneFeatures.internal(anno2)
  dot = crossprod(v1,v2)	
}


# calculate term similarities for a list of GO terms
getTermSim<-function(termlist, method="relevance", verbose=FALSE){
  S<-matrix(0,nrow=length(termlist),ncol=length(termlist))
  colnames(S)<-termlist
  rownames(S)<-termlist
  if(method %in% c("diffKernel")){		
    K = mget("K", envir=GOSimEnv, ifnotfound=list(function(x) stop(paste("Diffusion kernel not loaded!\nPlease invoke load.diffusion.kernel().", sep = ""), call. = FALSE)))$K
    Ktmp = matrix(NA, ncol=length(termlist), nrow=length(termlist))
    diag(Ktmp) = 1
    termlist = intersect(termlist, colnames(K))
    K = K[termlist, termlist]		
    return(K)
  }
  for(i in 1:length(termlist)){
    S[i,i] <- calcTermSim(termlist[i],termlist[i], method, verbose)				
    if(i > 1){
      for(j in 1:(i-1)){				
        S[i,j]<- calcTermSim(termlist[i],termlist[j], method, verbose)				
        S[j,i]<-S[i,j]
      }
    }
  }
  S
}


# get GraSM common disjunctive ancestors of two terms
getDisjCommAnc<-function(term1, term2){
  if(!exists("GOSimEnv")) initialize()
  ancestor<-get("ancestor",envir=GOSimEnv)
  IC<-get("IC", envir=GOSimEnv)
  if(term1 == term2){
    return(term1)
  }
  else{
    an1<-unlist(ancestor[names(ancestor) == term1])
    an2<-unlist(ancestor[names(ancestor) == term2])
    case1<-which(an2 == term1)  # term1 is an ancestor of term2
    case2<-which(an1 == term2) # term2 is an ancestor of term1
    if(length(case1) > 0){
      ancommon<-an1
      andisj<-getDisjAnc(term1, an1)
    }
    else if(length(case2) > 0){			
      ancommon<-an2
      andisj<-getDisjAnc(term2, an2)			
    }
    else{
      ancommon<-intersect(an1,an2)
      andisj<-getDisjAnc(c(term1,term2), ancommon) # we only need to calculate the disjunctives among the common ancestors!
      #andisj = unique(rbind(getDisjAnc(term1, an1), getDisjAnc(term2, an2)))
    }		
    djca<-c()
    cond1<-sapply(ancommon, function(x) setdiff(ancommon[which(IC[ancommon] >= IC[x])],x)) # which common ancestors are more informative than a_i?
    if(length(cond1) > 0){ # look for those that are disjunctive
      names(cond1)<-ancommon					
      for(i in 1:length(cond1)){
        res<-sapply(cond1[i][[1]],function(x) any(andisj[,1] == names(cond1[i]) & andisj[,2] == x) | any(andisj[,2] == names(cond1[i]) & andisj[,1] == x))				
        if(all(res))
          djca<-c(djca, names(cond1[i]))
      }
      djca= unique(djca)			
    }	
    if(length(djca)==0)
      djca<-ancommon[which.max(IC[ancommon])]# take minimum subsumer otherwise		
    return(djca)		
  }	
}

getGOGraph<-function(term, prune=Inf){
  #if(!require(graph))
  #	stop("Package graph is required for function getGOGraph")
  if(!exists("GOSimEnv")) initialize()	
  ontology<-get("ontology",envir=GOSimEnv)		
  if(ontology == "BP")
    G<-GOGraph(term,GOBPPARENTS)
  else if(ontology == "MF")
    G<-GOGraph(term,GOMFPARENTS)
  else if(ontology == "CC")
    G<-GOGraph(term,GOCCPARENTS)
  else
    stop(paste("ontology", ontology, "not known!"))			
  if(prune != Inf){
    dis = johnson.all.pairs.sp(G)		
    inc = unique(unlist(sapply(term, function(t) names(dis[t,])[dis[t,] < prune])))
    G = subGraph(nodes(G)[inc], G)
  }		
  G
}

# get FuSsiMeg depth factor
getDepthFactor<-function(term,G){	
  #if(!require(RBGL))
  #	stop("Package RBGL is required for function getDepthFactor")
  if(!exists("GOSimEnv")) initialize()	
  d<-sp.between(G,term,"all")[[1]]$length + 1  # start with depth = 1!
  D<-((d+1)/d)^get("alphaParam",envir=GOSimEnv)
  D
}

# get FuSSiMeg density factor
getDensityFactor<-function(term){
  if(!exists("GOSimEnv")) initialize()
  nchildren<-get("nchildren",envir=GOSimEnv)
  nparents<-get("nparents",envir=GOSimEnv)
  e<-nchildren[term] + nparents[term]	  
  betaParam<-get("betaParam",envir=GOSimEnv)
  E<-(1-betaParam)*get("Eavg",envir=GOSimEnv)/e + betaParam
  E
}

getGeneFeatures.internal = function(anno){
  ancestor<-get("ancestor",envir=GOSimEnv)	
  an<-unlist(ancestor[names(ancestor) %in% anno])	
  IC<-get("IC", envir=GOSimEnv)
  v = double(length(IC))	
  names(v) = names(IC)	
  v[c(anno, an)] = IC[c(anno, an)]
  v = v[!is.na(v)]
  v
}

# get GraSM disjunctive ancestors of a set of terms with ancestors an
getDisjAnc<-function(term, an){		
  #if(!require(RBGL))
  #	stop("Package RBGL is required for function getDisjAnc")
  G<-getGOGraph(term)	
  disan<-matrix(0,ncol=2,nrow=0)
  for(n1 in 1:length(an)){
    if(n1 > 1){
      for(n2 in 1:(n1-1)){
        if(!separates(term, an[n1], an[n2], G) && !separates(term, an[n2], an[n1], G))
          disan<-rbind(disan,c(an[n1], an[n2]))
      }
    }
  }
  disan
}

GOGraph = function(term, env){
  oldEdges <- vector("list", length = 0)
  oldNodes <- vector("character", length = 0)
  newN <- term
  done <- FALSE
  while (!done) {
    newN <- newN[!(newN %in% oldNodes)]
    if (length(newN) == 0)
      done <- TRUE
    else {
      oldNodes <- c(oldNodes, newN)
      numE <- length(newN)
      nedges <- AnnotationDbi::mget(newN, env = env, ifnotfound = NA)
      nedges <- nedges[!is.na(nedges)]
      oldEdges <- c(oldEdges, nedges)
      if (length(nedges) > 0)
        newN <- sort(unique(unlist(nedges)))
      else newN <- NULL
    }
  }
  rE <- vector("list", length = length(oldNodes))
  names(rE) <- oldNodes
  rE[names(oldEdges)] <- oldEdges
  rE <- lapply(rE, function(x) match(x, oldNodes))
  names(oldNodes) = oldNodes
  return(new("graphNEL", nodes = oldNodes, edgeL = lapply(rE,
              function(x) list(edges = x)), edgemode = "directed"))
}
