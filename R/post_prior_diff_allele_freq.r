post_prior_diff_allele_freq=function(sillyvec,loci,mydir,nchains,mixname){

  require("outliers")

  baseline=calc_freq_pop(sillyvec,loci)

  loci=dimnames(baseline[1,,])[[1]]

  nalleles=apply(baseline[1,,],1,function(vec){sum(!is.na(vec))})

  C=dim(baseline)[1]

  L=dim(baseline)[2]

  rawQ=baseline

  for(i in 1:C){
    for(loc in 1:L){
       rawQ[i,loc,1:nalleles[loc]]=(baseline[i,loc,1:nalleles[loc]]+1/nalleles[loc])/sum(baseline[i,loc,1:nalleles[loc]]+1/nalleles[loc]) 
    }
  }

  ChainNames=paste("Chain",1:nchains,sep="")

  mydirs=paste(mydir,"\\",mixname,ChainNames,"FRQ.FRQ",sep="")

  myfiles=lapply(mydirs,function(dr){myfile=readLines(dr);myfile[(length(myfile)/2+1):length(myfile)]})

  names(myfiles)=ChainNames

  myQ=lapply(myfiles,function(myfile){
        ans=t(vapply(myfile,function(strng){
              vec=strsplit(strng," ")[[1]];
              vec=vec[vec!=""];
              vec2=rep(NA,max(nalleles)+4);
              vec2[1:length(vec)]=vec;
              as.numeric(vec2)
            },as.numeric(rep(NA,4+max(nalleles)))));
        dimnames(ans)=list(1:nrow(ans),c("Iteration","Pop","Locus","n",paste("Allele",1:max(nalleles),sep="")));
        ans=data.frame(ans);
        q=rawQ;
        for(i in 1:C){
          for(loc in 1:L){
            q[i,loc,1:nalleles[loc]]=apply(ans[ans$Pop==i & ans$Locus==loc,5:(4+nalleles[loc])],2,mean)
          }
        };
        q;
  })

  names(myQ)=ChainNames

  myabsdiffsByPopByLocus=sapply(myQ,function(qq){ sapply(loci,function(locus){ apply(abs(qq[1:C,locus,1:nalleles[locus]]-rawQ[1:C,locus,1:nalleles[locus]]),1,mean) }) },simplify=FALSE)

  myabsdiffsByPop=sapply(myabsdiffsByPopByLocus,function(diffs){apply(diffs,1,mean)},simplify=FALSE)

  Outliers=sapply(myabsdiffsByPop,function(popdiffs){outlier(popdiffs,logical=TRUE)})

  return(list(Outliers=Outliers,MeanAbsDiffByPop=myabsdiffsByPop,PriorAlleleFreqs=rawQ,PosteriorAlleleFreqs=myQ))
}







