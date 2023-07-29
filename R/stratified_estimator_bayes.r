#' Stratified Estimator using Bayesian output
#'
#' This function reads in BAYES output files (all chains) from multiple temporal strata (mixtures) and, by weighting by harvest from each stratum, produces a seasonal total harvest 
#' estimate and C.I. for all defined regions (groups). Appologies for the number of arguments--this is necessary to give the required flexibility. The last four arguments will rarely need to be changed.
#'
#' WARNING -- A certain directory structure is required for successful execution. Particularly, each mixture's BAYES output files must be within their own sub-directory with the exact name as the mixture.
#' Also, the BAYES output files need names that begin with the mixture name, followed by the the prior name, followed by Chain -- e.g. "MixnamePriorChain1RGN.RGN" or "MixnamePriorChain1BOT.BOT".
#' The exact same prioname name must be used in all BAYES output file names.
#'
#' The output is an Excel file summarizing total harvest.
#'
#' @param groupvec A vector, the same length as the number of populations, each element is the group number for the corresponding population.
#'                e.g. groupvec = c(1,1,1,2,2,3) will combine the first 3 pops into group 1, the 4th & 5th pops into group 2 and the 6th pop into group 3.
#'                The elements of groupvec do not need to be in order -- e.g. groupvec = c(1,3,1,2,2,3) is acceptable.
#' @param groupnames Character string vector of names for the groups.
#' @param maindir The pathname to the sub-directories that contain the BAYES output files. This is also where the Excel output file will be written to.
#' @param mixvec Character string vector of mixture names, one for each stratum.
#' @param catchvec Vector of harvest numbers. Must have the same length as mixvec. Each element is the harvest for its corresponding stratum (mixture).
#' @param newname Chosen name of the output file without extension.
#' @param priorname Character string giving the name (type) of the prior used.
#' @param ext Used to indicate if "RGN" or "BOT" files are to be read in. Defaults to "RGN".
#' @param nchains The number of chains to read in from each mixture. Defaults to 3.
#' @param burn Proportion of chain length to use as burn-in. Defaults to 0.5.
#' @param alpha A proportion. Used to create a 100*(1-alpha)% confidence interval. Defaults to 0.1 (90% CI).
#' @param xlxs Logical; if TRUE, the output will be written to an Excel file. Requires the `xlsx` package.
#' @param PosteriorOutput Logical; if TRUE, the function will return both the summary statistics and the posterior output. If FALSE, it will return only the summary statistics.
#'
#' @return A summary of the stratified estimator, including mean, standard deviation, median, 90% confidence interval, proportion of values less than the threshold, and Gelman-Rubin statistic.
#' @examples
#' \dontrun{
#' my.groupvec <- c(1,1,1,1,1,1,1,1,2,3,5,5,4,5) 
#' my.groupnames <- c("Alaska","Nass","Skeena","Fraser","Other") 
#' my.maindir <- "V:\\WORK\\SARA\\Sockeye\\SEAK\\Analysis\\Mixtures\\0405Dists101104\\Bayes\\Feb2009\\Output\\2005dist104" 
#' my.mixvec <- c("year2005dist104stat2829a","year2005dist104stat30a","year2005dist104stat31a","year2005dist104stat32a","year2005dist104stat33a","year2005dist104stat34a") 
#' my.catchvec <- c(14904,20786,36201,93584,66464,274455)
#' my.newname <-  "TotalHarvest2005dist104_90percentCI"
#' my.priorname <- "Flat" 
#'
#' StratifiedEstimator <- stratified_estimator_bayes(groupvec=my.groupvec, groupnames=my.groupnames, maindir=my.maindir, mixvec=my.mixvec, catchvec=my.catchvec, newname=my.newname, priorname=my.priorname)
#' }
#'
#' @export
#' 
stratified_estimator_bayes=function(groupvec, groupnames, maindir, mixvec, catchvec, CVvec=rep(0,length(catchvec)), newname, priorname="", ext="RGN", nchains=5, burn=0.5, alpha=0.1, threshold=5e-7, xlxs=TRUE, PosteriorOutput=FALSE){

nstrata=length(mixvec)

if(nstrata==1){stop("This function is for multiple strata!!!")}

if(length(mixvec) != length(catchvec)) {stop("mixvec and catchvec must be the same length")}
   
if(xlxs){
  while(!require(xlsx)){install.packages("xlsx")}
}

while(!require(coda)){install.packages("coda")}

chains=paste("Chain",1:nchains,sep="")

files=sapply(mixvec,function(mix){sapply(chains,function(chain){paste(maindir,"\\",mix,"\\",mix,priorname,chain,ext,".",ext,sep="")},simplify=FALSE)},simplify=FALSE)

output00=rapply(files,read.table,row.names=1,how = "list")

nsamps00=BurnIn=nsamps0=nsamps=array(NA,c(nchains,nstrata),list(chains,mixvec))
 
output0=Keep=output=setNames(vector("list",nchains),chains)

for(chain in chains){
  
  output0[[chain]]=setNames(vector("list",nstrata),mixvec)

  for(mix in seq_along(mixvec)){
    
    nsamps00[chain,mix]=nrow(output00[[mix]][[chain]])

    BurnIn[chain,mix]=floor(burn*nsamps00[chain,mix])

    output0[[chain]][[mix]]=output00[[mix]][[chain]][-(1:BurnIn[chain,mix]),]

    nsamps0[chain,mix]=nrow(output0[[chain]][[mix]])

  }

}


Thin=floor(nsamps0/min(nsamps0))

for(chain in chains){
  
  output[[chain]]=setNames(vector("list",nstrata),mixvec)

  Keep[[chain]]=setNames(vector("list",nstrata),mixvec)

  for(mix in seq_along(mixvec)){
    
    ans0=1:nsamps0[chain,mix]

    Keep[[chain]][[mix]]=ans0[ans0%%Thin[chain,mix]==0]

    nsamps[chain,mix]=length(Keep[[chain]][[mix]])

    output[[chain]][[mix]]=output0[[chain]][[mix]][Keep[[chain]][[mix]],]

  }

}

nCombinedSamps=apply(nsamps,2,sum)

if(sd(nCombinedSamps)){

  message("Not all mixtures were run for same number of iterations - subsampling to mininum number of iterations in mixvec:")
  
  prob_mix <- colnames(nsamps)[apply(nsamps, 2, function(x) { min( x ) > min( nsamps )})]
  
  message(paste(prob_mix))
  
  set.seed(12345) # seed for reproducability
  subsampler <-
    sample(
      x = 1:max(nsamps),
      size = min(nsamps),
      replace = FALSE
    ) # subsampler chooses values (x) from largest dataframe, then downsamples to smallest.
  for (chain in chains) {
    for (mix in prob_mix) {
      output[[chain]][[mix]] = output[[chain]][[mix]][subsampler,]
      nsamps[chain, mix] = nrow(output[[chain]][[mix]])
    }
  }
  nCombinedSamps = apply(X = nsamps, MARGIN =  2, FUN =  sum)
}

n=unique(nCombinedSamps)

groupvec=as.numeric(groupvec)
catchvec=setNames(as.numeric(catchvec),mixvec)
CVvec=setNames(as.numeric(CVvec),mixvec)

C=length(groupvec)

G=max(groupvec)

WghtAve0=WghtAve=H0=setNames(vector("list",nchains),chains)

for(chain in chains){

  H0[[chain]]=sapply(seq_along(mixvec),function(mix){lnvar=log(CVvec[mix]^2+1);lnmean=log(catchvec[mix])-lnvar/2;rlnorm(nsamps[chain,mix],lnmean,sqrt(lnvar))},simplify=FALSE)

  WghtAve0[[chain]]=Reduce("+",lapply(seq_along(mixvec),function(mix){H0[[chain]][[mix]]*output[[chain]][[mix]]}))/Reduce("+",lapply(seq_along(mixvec),function(mix){H0[[chain]][[mix]]}))

  WghtAve[[chain]]=t(rowsum(t(WghtAve0[[chain]]),group=groupvec))

  colnames(WghtAve[[chain]]) <- groupnames
  
  WghtAve[[chain]]=as.mcmc(WghtAve[[chain]])

}

WghtAve=as.mcmc.list(WghtAve)

OutputWghtAve=Reduce(rbind,WghtAve)

summary=array(NA,c(G,7),list(groupnames,c("mean", "sd", "median", paste0(100*alpha/2,"%"), paste0(100*(1-alpha/2),"%"), "P=0", "GR")))

summary[groupnames,"mean"]=apply(OutputWghtAve,2,mean)

summary[groupnames,"sd"]=apply(OutputWghtAve,2,sd)

summary[groupnames,"median"]=apply(OutputWghtAve,2,median)

summary[groupnames,paste0(100*alpha/2,"%")]=apply(OutputWghtAve,2,quantile,probs=alpha/2)

summary[groupnames,paste0(100*(1-alpha/2),"%")]=apply(OutputWghtAve,2,quantile,probs=1-alpha/2)

summary[groupnames,"P=0"]=apply(OutputWghtAve,2,function(clm){sum(clm<threshold)/n})

summary[groupnames,"GR"]=gelman.diag(WghtAve,multivariate=FALSE,transform=TRUE)[[1]][,1]

if(xlxs){
  write.xlsx(summary, file=paste(maindir, "\\", newname, ".xlsx", sep=""))
}

if(PosteriorOutput) {
  ans = list(Stats=summary,Output=OutputWghtAve)
} else {
  ans = summary
}

return(ans)

}







