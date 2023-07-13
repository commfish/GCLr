read_biomark_qc=function(qccsvFilepaths,skip=15){
###################################################################################################################################################################################################################################################
#
#qccsvFilepaths=c("V:/DATA/All SNP data/Chinook/Project K42 Chinook Inseason_Postseason 2012/Project K42 M192 BioMark/qc/Data/Single Chip/K42 Combined Output.csv","V:/DATA/All SNP data/Chinook/Project K45- AYK, BB Chinook baseline/Project K45 192chip BioMark/qc/Data/Single Chip/K45_qc2_1381716192/K45_qc2_ModifiedforR.csv")
#
#
#
###################################################################################################################################################################################################################################################

  sillyvec=as.vector(sapply(objects(pattern="*\\.gcl",pos=1),function(gclname){strsplit(gclname,split="\\.gcl")[[1]][1]}))

  Genotypesqc=NULL

  for(pth in qccsvFilepaths){
    Genotypesqc=rbind(Genotypesqc,read.table(pth,header=TRUE,sep=",",colClasses="character",stringsAsFactors=FALSE,skip=skip)[,c("Name","Converted","Assay","ID","Allele.X","Allele.Y")])
  }
  
  Genotypesqc=Genotypesqc[Genotypesqc$Name!="NTC",] #Added the removal of NTCs 10/2/2012 A.B.
  
  Genotypesqc$Converted[Genotypesqc$Converted=="No Call"]=NA
  
  Genotypesqc$Converted[Genotypesqc$Converted=="Invalid"]=NA

  MyTab=table(Genotypesqc$Name,Genotypesqc$Assay)

  SillyVialqc=dimnames(MyTab)[[1]]

  nVials=length(SillyVialqc)
  
  lociqc=dimnames(MyTab)[[2]]

  nlociqc=length(lociqc)

  loci=LocusControl$locusnames

  nalleles=LocusControl$nalleles[loci]

  ploidy=LocusControl$ploidy[loci]

  nloci=length(loci)

  alleles=LocusControl$alleles[loci]

  names(loci)=loci

  IND=match(lociqc,loci)

  lociqcmissng=lociqc[is.na(IND)]

  if(sum(is.na(IND))){

    print(paste("qc loci: ",paste("'",lociqcmissng,"'",collapse=",")," do not exsist in 'LocusControl'!!! : Carrying On...!!!",sep=""))

    lociqc=lociqc[!is.na(IND)]

    nlociqc=length(lociqc)
    
  }


  Scores=array(NA,c(nVials,nloci),dimnames=list(SillyVialqc,loci))

  for(vial in SillyVialqc){
    
    for(locus in loci){
     
      Scores[vial,locus]=Genotypesqc[Genotypesqc$Name==vial & Genotypesqc$Assay==locus,"Converted"][1]
   
    }
  
  }
  

  Genotypesqcsillys=as.character(sapply(SillyVialqc,function(nm){strsplit(nm,split="_")[[1]][1]}))

  sillyvecqc0=unique(Genotypesqcsillys)

  sillyvecqc=sillyvecqc0[!is.na(match(sillyvecqc0,sillyvec))]

  GenotypesqcIDs=as.character(sapply(SillyVialqc,function(nm){strsplit(nm,split="_")[[1]][2]}))

  AttributsColNames=dimnames(get(objects(pattern="*\\.gcl",pos=1)[1],pos=1)$attributes)[[2]]

  for(silly in sillyvecqc){
    IND=Genotypesqcsillys==silly
    IDs=GenotypesqcIDs[IND]
    n=length(IDs)
    scores=array(NA,c(n,nloci,max(ploidy)),dimnames=list(IDs,loci,paste("Dose",1:max(ploidy),sep="")))
    counts=array(NA,c(n,nloci,max(nalleles)),dimnames=list(IDs,loci,paste("Allele",1:max(nalleles),sep="")))
    for(id in IDs){
      idsillyvial=paste(silly,id,sep="_")
      for(locus in loci){
        if(!is.na(Scores[idsillyvial,locus])){
          scores[id,locus,1:ploidy[locus]]=sort(strsplit(Scores[idsillyvial,locus],split=":")[[1]][1:ploidy[locus]])
          for(alll in 1:nalleles[locus]){
            counts[id,locus,alll]=sum(scores[id,locus,1:ploidy[locus]]==alleles[[locus]][alll])
          }
        }
      }
    }
    attributes=data.frame(array(NA,c(n,length(AttributsColNames)),dimnames=list(IDs,AttributsColNames)))
    attributes$FK_FISH_ID=IDs
    attributes$SillySource=as.character(sapply(SillyVialqc[IND],function(sllyvl){paste(strsplit(sllyvl,split="_")[[1]],collapse="qc_")}))
    assign(paste(silly,"qc.gcl",sep=""),list(counts=counts,scores=scores,n=n,attributes=attributes),pos=1)
  }
assign(x = "qcSillys", value = paste(sillyvecqc,"qc",sep=""), pos = 1)
# return(paste(sillyvecqc,"qc",sep=""))
}






                         
