var_comp <- function(sillyvec, loci, groupvec, fstatfile = NULL, dir){
  ####################################################################################################################################################################################################################################################################
  #
  # This function provides the variance componets of a 3 level hierarchy
  # See GDA Weir, 1995. pg 184
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # sillyvec = Kodiak49Pops
  #   ~ The above definitions in () assume that you are running this on sillys that are pops
  #   ~ This could also be done on sillys that are collections, but the meaning would be different
  # loci = loci90
  #   ~ Can include haploid markers and/or combined loci
  # groupvec = Kodiak49GroupVec9
  # file = "fstatfile.dat"
  #   ~ If you already have an fstat.dat file, specify here, otherwise it defaults to NULL
  # dir = "V:/WORK/Sockeye/Kodiak/Kodiak Afognak Baseline Look 2014/Kyle/ThetaP"
  #   ~ This is the location that the function will first search for an "fstatfile.dat", if none is found, it will create one and put it there
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function is intended to help with marker selection (i.e. going from a 96 -> 24 SNP baseline)
  # Briefly the four variance components in the output (matrix) are:
  #   P - for populations (RGs in our case)
  #   S - for subpopulations within populations (populations within RGs in our case)
  #   I - for individuals within subpopulations (individuals within populations in our case)
  #   G - for alleles within individuals (akin to heterzygosity)
  # 
  #
  # Created by Kyle Shedd on Thu Jul 30 15:26:28 2015
  # Original function is "BackwardEliminationThetaP.GCL" in the temp folder
  #
  ####################################################################################################################################################################################################################################################################
    
  require("hierfstat")
  
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