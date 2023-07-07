#' @title Create HWLER Control Files
#'
#' @description This function creates HWLER control files for each chain based on the specified parameters.
#'
#' @param sillyvec A character vector of population sillys used to create the baseline file.
#' @param loci  character vector of the loci used to produce the baseline and mixture files.
#' @param input The base name for the control files.
#' @param mixbase A character vector indicating mixture and/or baseline individuals are present in the data.
#' @param nsamples A vector specifying the number of samples for each chain.
#' @param nchains The number of MCMC chains to analyze the mixtures.
#' @param dir The directory path where the control files will be saved.
#' @param initval The initial starting value for each MCMC chain.
#' @param seeds A matrix of random seeds containing 3 seeds per chain, where nrow(seeds) == 3 and ncol(seeds) == nchains.
#' @param thin Thinning intervals for MCMC sample of 1) stock proportions, 2) baseline allele or type relative frequencies, and 3) stock assignments of each mixture individual.
#' @param inputfortran The input is in The Fortran format.
#' @param switches A character string of logical switches, with default value "T T T F T T T T T" (see details).
#' 
#' @details The `switches` argument has 9 program options turned “on” with “T” for true or “off” with “F” for false:
#'   \enumerate{
#'     \item Print summary of baseline data to summary file
#'     \item Print mixture data to summary file
#'     \item Output MCMC population assignments of individuals or snapshots of cluster assignments of individuals 
#'     \item Suppress the output of MCMC samples of population proportions or snapshots of numbers of individuals assigned to each cluster 
#'     \item Saseline samples provided
#'     \item Mixture samples provided
#'     \item Alphanumeric identification for each baseline or mixture individual provided at the end of each line of data file
#'     \item Output MCMC samples of partitions for making a binary tree
#'     \item Perform MCMC sampling of the Dirichlet mass parameter
#'     }
#'     
#' @return Writes out HWLER control (.ctl) files.
#'
#' @examples
#' GCLr::create_hwler_ctl()
#'
#' @export
create_hwler_ctl<-function(sillyvec,loci,input,mixbase="mix",nsamples=c(1000,2,5),
                               nchains=5,dir,initval="1.0",seeds=matrix(sample(seq(10000),3*nchains),
                               nrow=3),thin=c(1,1,1),inputfortran,switches="T T T F T T T T T"){


  chains=paste("Chain",1:nchains,sep="")

  dirs=paste(dir,"\\",input,"Chain",1:nchains,".ctl",sep="")
  names(dirs)=chains

  nsillys=length(sillyvec)

  nalleles=LocusControl$nalleles[loci]

  nloci=length(loci)
  
  seeds=cbind(seeds)
  dimnames(seeds)=list(1:3,chains)

  files=lapply(chains,function(chain){paste(input,chain,sep="")})
  names(files)=chains

  files=lapply(files,function(file){rbind(file,paste(input,mixbase,sep="."))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(input,chain,"SUM.SUM",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(input,chain,"BOT.BOT",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(input,chain,"CLS.CLS",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(input,chain,"STK.STK",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(input,chain,"TXT.TXT",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(input,chain,"ALP.ALP",sep=""))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],paste(nsamples[1],nsamples[2],nsamples[3],colapse="    "))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],length(sillyvec))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],length(loci))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],initval)})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],cbind(seeds[,chain]))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],cbind(thin))})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],inputfortran)})
  names(files)=chains

  files=lapply(chains,function(chain){rbind(files[[chain]],switches)})
  names(files)=chains


  files=lapply(seq(length(chains)),function(chain){rbind(files[[chain]],cbind(sapply(1:nloci,function(d){paste(sprintf("%3s",d),sprintf("%2s",nalleles[loci[d]]),"",loci[d],collapse="")})))})
  names(files)=chains


  
  files=lapply(chains,function(chain){rbind(files[[chain]],cbind(sapply(1:nsillys,function(i){paste(sprintf("%3s",i),"", format(ifelse(nchar(sillyvec[i])<18,sillyvec[i],substr(sillyvec[i],start=1,stop=18)),width=18),collapse="")})))})
  names(files)=chains

  empty=sapply(chains,function(chain){write.table(files[[chain]],dirs[chain],quote=FALSE,row.names=FALSE,col.names=FALSE)})

  return(NULL)


}
