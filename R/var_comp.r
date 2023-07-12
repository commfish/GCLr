#' Variance Components Calculation
#'
#' This function calculates the variance components of a three-level hierarchy. This function is intended to help with marker selection (i.e., going from a 96 to 24 SNP baseline).
#' 
#' @param sillyvec A character vector specifying the population silly codes for analysis.
#' @param loci A character vector of locus names to include in the analysis. Can include haploid markers and/or combined loci.
#' @param groupvec A numeric vector specifying the reporting group affiliation for each silly code.
#' @param fstatfile Optional character string specifying the path to an existing FSTAT file (".dat"). If provided, the function verifies compatibility with the input data. Default is NULL.
#' @param dir A character string specifying the directory where an FSTAT file (".dat") will be created. This argument is used only when `fstatfile` is NULL.
#' 
#' @details 
#' The original function is "BackwardEliminationThetaP.GCL" in the temp folder.
#' 
#' @return A matrix containing the four variance components for each locus:
#'   \itemize{
#'     \item P: For populations (reporting groups in this case).
#'     \item S: For subpopulations within populations (populations within reporting groups).
#'     \item I: For individuals within subpopulations (individuals within populations).
#'     \item G: For alleles within individuals (akin to heterozygosity).
#'    }
#'  Component definitions in () assume that you are running this function on silly codes that are populations. This function could also be run on silly codes that are collections, but the definitions in () would be different.
#' 
#' @examples
#' GCLr::var_comp(sillyvec = Kodiak49Pops, loci = loci90, groupvec = Kodiak49GroupVec9, dir = "V:/WORK/Sockeye/Kodiak/Kodiak Afognak Baseline Look 2014/Kyle/ThetaP")
#' 
#' @references Weir, B.S. 1995. Genetic Data Analysis. Sinauer Associates, Inc. Sunderland, MA. (pg 184)
#' 
#' @export
var_comp <- function(sillyvec, loci, groupvec, fstatfile = NULL, dir){

  names(sillyvec) <- NULL
  
  nloci=length(loci)
  nalleles=LocusControl$nalleles[loci]
  ploidy=LocusControl$ploidy[loci] 
  alleles=LocusControl$alleles[loci]
  
  my.gcl=sapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)},simplify=FALSE)
  
  n=sapply(sillyvec,function(silly){my.gcl[[silly]]$n})
  
  if(!is.null(fstatfile)) {
    
    dat0=read.fstat.data(fstatfile)
    
    if(dim(dat0)[1] != sum(n)) {stop(paste(fstatfile, " does not have the same number of fish as 'sillyvec', verify that the 'fstatfile' corresponds to the fish you want to analyze!"))}

    if(dim(dat0)[2] != (nloci + 1)) {stop(paste(fstatfile, " does not have the same number of loci as 'loci', verify that the 'fstatfile' corresponds to the loci you want to analyze!"))}
    
    if(length(unique(dat0$Pop)) != length(groupvec)) {stop("Number of pops in 'fstatfile' is not equal to length(groupvec)!")}
    
  } else {
    cat("No 'fstatfile' provided, so creating one\n")
    
    fstatdir=paste(dir,"\\", length(groupvec), "pops", nloci, "loci", "_","fstatfile.dat",sep="")
        
    maxchar=nchar(nalleles)+1
    names(maxchar)=loci
    
    nsillys=length(sillyvec)
    
    maxsillychar=nchar(nsillys)+1
    
    scores=sapply(sillyvec,function(silly){my.gcl[[silly]]$scores[,loci,]},simplify=FALSE)
    
    counts=sapply(sillyvec,function(silly){
      sapply(1:n[silly],function(i){
        paste(c(match(silly,sillyvec),sapply(loci,function(locus){
          ifelse(is.na(scores[[silly]][i,locus,1]),paste(rep(0,ploidy[locus]*maxchar[locus]),collapse=""),paste(sapply(1:ploidy[locus],function(allele){
            paste(c(rep(0,maxchar[locus]-nchar(match(scores[[silly]][i,locus,allele],alleles[[locus]]))),match(scores[[silly]][i,locus,allele],alleles[[locus]])),collapse="")
          }),collapse=""))
        })),collapse=" ")
      })
    },simplify=FALSE)      
    
    fstat=paste(nsillys,nloci,max(nalleles),max(maxchar),sep=" ")
    
    fstat=rbind(fstat,cbind(loci))
    
    fstat=rbind(fstat,cbind(as.vector(unlist(counts))))
    
    write.table(x = fstat, file = fstatdir, row.names = FALSE, col.names = FALSE, quote = FALSE) 
    
    dat0=read.fstat.data(fstatdir)
    
  }
  
  dat=data.frame(rep(groupvec,n),dat0)
  
  dimnames(dat)[[2]]=c("Region","Pop",loci)
  
  G=max(groupvec)
  
  VarCompByLocus=array(0,c(nloci,4),dimnames=list(loci,c("P","S","I","G")))
  
  cat("\nCalculating Variance Components for each locus\n", sep = '')
  if (.Platform$OS.type == "windows") flush.console()
  pb <- txtProgressBar(min = 0, max = nloci, style = 3)
  
  for(locus in loci){
    setTxtProgressBar(pb = pb, value = which(loci == locus))
    
    data=data.frame(dat[,1:2],dat[,locus])
    
    vc=varcomp(data,diploid=ploidy[locus]>1)$overall
    
    if(!sum(is.na(vc))){
      
      VarCompByLocus[locus,1:length(vc)]=vc
      
    }
    
  }
  
  return(VarCompByLocus)
}