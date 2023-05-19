custom_comb_hwler_output=function(groupvec, groupnames, maindir, mixvec, ext="BOT", burn=0.5, alpha=0.1,PosteriorOutput=TRUE){
#######################################################################################################################################################################################################################################################################################################################################################################################################################################
#
#    Notice:  This function makes all accomodations for the extra baseline group.
#
#    "groupvec" is the same length as the number of populations in the KNOWN baseline.
#
#    "groupnames" is the same length as the number of groups in the KNOWN baseline.
#  
#
#  
#
# 
#
#  
#
#######################################################################################################################################################################################################################################################################################################################################################################################################################################

  G=max(groupvec)

  groupvec=c(groupvec,G+1)

  groupnames=c(groupnames,"GroupX")

  Results=vector("list",length(mixvec))
  names(Results)=mixvec

  pdf(file=paste(maindir,"\\TracePlotsHWLER",format(Sys.time(),"%b%d%Y_%H_%M_%S"),".pdf",sep=""),width=11, height=8, family="Helvetica",pointsize=20)

  for(mix in mixvec){

    HWLERoutput=read.table(paste(maindir,"\\",mix,"Chain1",ext,".",ext,sep=""))[,-1]

    nits=nrow(HWLERoutput)

    Burn=max(1,floor(burn*nits))

    GroupedHWLERoutput=array(NA,c(nits,G+1),dimnames=list(1:nits,groupnames))

    for(group in groupnames){

      g=match(group,groupnames)

      if(sum(groupvec==g)>1){

        GroupedHWLERoutput[1:nits,group]=apply(HWLERoutput[1:nits,groupvec==g],1,sum)
     
      }      

      if(sum(groupvec==g)==1){

        GroupedHWLERoutput[1:nits,group]=HWLERoutput[1:nits,groupvec==g]
     
      } 
   
      plot(GroupedHWLERoutput[1:nits,group],type="l",ylim=c(0,1),yaxp=c(0,1,10),col="red",xlab="Iteration",ylab="Proportion",main=paste("Trace of ",group," from ",mix,sep=""))

      abline(v=Burn,col="blue")     

    }

    UsedHWLERoutput=GroupedHWLERoutput[(Burn+1):nits,groupnames]

    Results[[mix]]=vector("list",PosteriorOutput+1)

    if(PosteriorOutput){
   
      names(Results[[mix]])=c("SummaryStats","PosteriorOutput")

      Results[[mix]][["SummaryStats"]]=data.frame(mean=apply(UsedHWLERoutput,2,mean),sd=apply(UsedHWLERoutput,2,sd),LowCI=apply(UsedHWLERoutput,2,function(clmn){quantile(clmn,alpha/2)}),UpCI=apply(UsedHWLERoutput,2,function(clmn){quantile(clmn,1-alpha/2)})) 
   
      Results[[mix]][["PosteriorOutput"]]=UsedHWLERoutput
    } 

    if(!PosteriorOutput){
   
      names(Results[[mix]])=c("SummaryStats")

      Results[[mix]][["SummaryStats"]]=data.frame(mean=apply(UsedHWLERoutput,2,mean),sd=apply(UsedHWLERoutput,2,sd),LowCI=apply(UsedHWLERoutput,2,function(clmn){quantile(clmn,alpha/2)}),UpCI=apply(UsedHWLERoutput,2,function(clmn){quantile(clmn,1-alpha/2)})) 
   
    }

    
   
  }

  dev.off()

  return(Results) 

}







