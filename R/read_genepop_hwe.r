read_genepop_hwe=function(file, sillyvec = NULL, summaryAsNumeric = TRUE){
  #############################################################################################################################
  #This function reads in the population output from a GENEPOP Hardy-Weinberg test ("*.P") file and returns a list of 2:
  #1) is a data frame containing all of the HWE data for each locus by Pop
  #2) is a matrix containing just the PValues for each Pop and the overall PValue across Pops (Fisher's method)
  #
  # "file" - the full file path, including the ".P" extension.  Make sure the file has not been modified.
  # "sillyvec" - optional character vector to provide for pop names as opposed to getting form genepop .P file
  # "summaryAsNumeric" - logical, indicates if you want summary p-values from exact test to be numeric or character
  #   Note: forcing to numeric takes output such as "High. sign." and forces to 0.000 for Genepop 4.6 and below
  #   Note: forcing to numeric takes output such as " < 0.1555" and forces to 0.1555 for Genepop 4.7 and above
  #
  # Example: HWE=read_genepop_hwe(file="V:/WORK/Chum/AHRP/Data/AHRPsamples.txt.P")
  #
  # Written 8/6/14 Kyle Shedd - data.frame with p-values for each pop + overall
  # Updated 8/22/14 08:16 to handle exact test vs. MCMC method
  # Updated 9/16/14 10:25 KS changed how DF are determined for Fisher's method (global p-values for loci across populations)
  # Updated 9/16/14 11:33 KS Global p-values for loci across populations (summary[,"Overall Pops"]) are now read directly from the 
  #  Genepop *.P file rather than calculated with Fischer's method in R. The reasoning for this new approach is that 1), it no 
  #  longer matters if you use the Exact test or MCMC method, and 2) since Genepop rounds p-values to the nearest 0.0001, small 
  #  p-values were simply displayed as 0.0000, which made subsequent calculations of global p-values via Fischer's method fail, 
  #  since you can't take the ln(0).
  # Updated 9/17/14 11:51 KS Added global p-values for populations across loci (summary ["Overall Loci",])
  # Updated 2/8/15 17:32 KS Added global p-value for Overall Pops over all loci (Fisher's Method). NOTE: this should only be
  #  trusted/used if the sample size per pop is similar...
  # Updated 7/27/15 20:00 KS Added an "if" "else" section to handle the calculation of "Overall Loci" and "Overall Pops"
  #  depending on whether p-values were calculated in Genepop with the "Exact" or "MCMC" setting. If MCMC then hard 0s are
  #  replaced with repzero = batches * iterations
  #############################################################################################################################
  
  while(!require(reshape2)){ install.packages("reshape2") }

  hwp=scan(file,what='',sep = "\n")

  npops=as.numeric(strsplit(hwp[grep("Number of populations detected:    ",hwp)],split="Number of populations detected:    ")[[1]][2])
  
  nloci=as.numeric(strsplit(hwp[grep("Number of loci detected:           ",hwp)],split="Number of loci detected:           ")[[1]][2])
  
  popstart=grep("Pop : ",hwp)
  
  if(is.null(sillyvec)) {
    
    pops=strsplit(hwp[popstart],split="Pop : ")
    
    pops=unlist(pops)[seq(2,2*npops,by=2)]
    
  } else {
    
    pops <- sillyvec
    
  }
  
  nmonoloci=length(grep("not diploid.",hwp))
  
  ndiploci=nloci-nmonoloci
  
  HWdata=NULL
  
  for(i in 1:npops){  
    HWdata=rbind(HWdata,cbind(Pop=rep(pops[i],ndiploci),colsplit(
      sapply((popstart[i]+6):(popstart[i]+6+ndiploci-1),function(row){gsub(pattern="[[:blank:]]+",x=hwp[row],replacement="/",fixed=F)}),
      pattern="/",names=c("Locus","PValue","SE","WC Fis","RH Fis","Steps","Switches","Low"))[,1:6]))
  }
  
  # Break if Pop = 1
  
  if(npops == 1) {print(cat("Only 1 POP detected, so no summary SummaryPValues.\nAll p-values are pulled directly from the Genepop *.P file, none are derived in this function!\nHaploid loci are not included in the output data.frame."))
                  return(HWdata)
  }
  
  # Proceed as normal
  
  loci=gsub(colsplit(hwp[grep("Locus ",hwp)],pattern="\"",names=c("Locus","Locus"))[,2],pattern="\"",replacement="")
  
  invisible(ifelse(length(grep(" not diploid.",loci))==0,assign("mono",0),assign("mono",grep(" not diploid.",loci))))
  
  invisible(ifelse(mono==0,assign("diploci",loci),assign("diploci",loci[-mono])))
  
  HWdata[,"Locus"]=diploci
  
  HWdata$PValue=suppressWarnings(as.numeric(gsub(HWdata$PValue,pattern="No",replacement="NA")))
  
  HWdata$SE=suppressWarnings(as.numeric(gsub(HWdata$SE,pattern="information.",replacement="NA")))
  
  loci=gsub(loci,pattern=" not diploid.",replacement="")
  
  smmry=matrix(nrow=nloci+1,ncol=npops+1,dimnames=list(c(loci,"Overall Loci"),c(pops,"Overall Pops")))
  
  for(i in 1:npops){
    for(j in 1:nloci){
      smmry[j,i]=HWdata[(match(x=loci[j],table=HWdata[,"Locus"])+((i-1)*ndiploci)),"PValue"]      
    }
  }
  
  ## Overall probability using Fisher's method, checked with Genepop to verify results, not the same when p-values = 0 for a population at a locus (rounding issue)
  # smmry[,"Overall"]=t(t(round(pchisq(-2*apply(log(smmry[,1:npops]),1,function(x) sum(x,na.rm=TRUE)),2*apply(smmry[,1:npops],1,function(x) sum(!is.na(x))),lower.tail=FALSE),4)))

  ## Overall probability using Fisher's method
  # If MCMC, replace hard 0 with repzero and calculate overall pops and overall loci via Fisher's ChiSqaure method
  # Else, if Exact test, p-values overall pops or overall loci are pulled directly from Genepop *.P file if Exact test
  
  if(length(grep(pattern = "Batches:", x = hwp)) == 1) {
    
    batches <- as.numeric(strsplit(hwp[grep("Batches:", hwp)], split = "Batches:                     ")[[1]][2])
    iterations <- as.numeric(strsplit(hwp[grep("Iterations per batch:        ", hwp)], split = "Iterations per batch:        ")[[1]][2])
    repzero <- as.numeric(1 / (batches * iterations))
    
    smmry[which(smmry == 0, arr.ind = TRUE)] = repzero
    
    overallpops <- round(pchisq(-2 * apply(log(smmry[1:nloci, 1:npops]), 1, function(locus) {sum(locus, na.rm = TRUE)} ), 2 * apply(smmry[1:nloci, 1:npops], 1, function(locus) {sum(!is.na(locus))} ), lower.tail = FALSE), 4)
    
    overallloci <- round(pchisq(-2 * apply(log(smmry[1:nloci, 1:npops]), 2, function(pop) {sum(pop, na.rm = TRUE)} ), 2 * apply(smmry[1:nloci, 1:npops], 2, function(pop) {sum(!is.na(pop))} ), lower.tail = FALSE), 4)
    
  } else {
    
    overallpops=colsplit(hwp[grep("Locus ",hwp)+5+npops+4],pattern=" Prob :    ",names=c("Trash","P-value"))[,2]
    
    overallloci=colsplit(hwp[grep("Pop : ",hwp)+5+ndiploci+4],pattern=" Prob :    ",names=c("Trash","P-value"))[,2]
    
    ## Force SummaryPValues matrix to be numeric for Exact test
    # If using Genepop version 4.6 or lower, "High. sign." is forced to "0.0000"
    # If using Genepop version 4.7 or higher, " < 0.0511" is forced to "0.0511"
    
    if(summaryAsNumeric) {
      
      if(any(grepl(pattern = "High. sign.", overallpops)) | any(grepl(pattern = "High. sign.", overallloci))) {
        
        overallpops=as.numeric(gsub(overallpops,pattern="High. sign.",replacement="0.0000"))
        
        overallloci=as.numeric(gsub(overallloci,pattern="High. sign.",replacement="0.0000"))
        
      }  # 4.6 and older
      
      if(any(grepl(pattern = "<", overallpops)) | any(grepl(pattern = "<", overallloci))) {
        
        overallpops=as.numeric(gsub(overallpops,pattern=" < ",replacement=""))
        
        overallloci=as.numeric(gsub(overallloci,pattern=" < ",replacement=""))
        
      }  # 4.7 and newer
      
    }  # summaryAsNumeric
    
  }  # Exact test
   
  smmry[1:nloci,"Overall Pops"]=overallpops
  
  smmry["Overall Loci",1:npops]=overallloci
  
  if(storage.mode(smmry) == "character") {
    
    smmry[nrow(smmry),ncol(smmry)] = NA
    
  } else {
    
    smmry[nrow(smmry),ncol(smmry)]=round(pchisq(q=(-2*sum(log(smmry[1:(nrow(smmry)-1),ncol(smmry)]))), df=(2*(nrow(smmry)-1)), lower.tail=FALSE), 4)
    
  }
    
  lst=list("DataByPop"=HWdata,"SummaryPValues"=smmry)
  
  print(cat("All p-values are pulled directly from the Genepop *.P file if calculated via Exact test.\nIf calculated via MCMC then over all loci/pops are derived in this function, correcting p = 0 for the number of batches * iterations!\nNA in 'smmry' dataframe means that locus in either 1) monomorphic (or very low MAF) or 2) not diploid (haploid mtSNP).\nHaploid loci are not included in the 'HWdata' data.frame."))
  return(invisible(lst))
}