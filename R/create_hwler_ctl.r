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
