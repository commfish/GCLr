#' @title Read Genepop HWE Output
#' 
#' @description
#' This function reads in output from a `genepop` Hardy-Weinberg test ("*.P") file.
#'
#' @param file The full file path to the `genepop` HWE output, including the ".P" extension.
#' @param sillyvec An optional character vector of silly codes without the ".gcl" extension (default = `NULL`).
#' If `NULL`, `Pop` names will come directly from the ("*.P") file and likely include "_fishID" extensions.
#' If supplying `sillyvec`, make sure it is the same `sillyvec` used in [GCLr::gcl2genepop()].
#' @param summaryAsNumeric A logical vector of length 1 indicating whether character p-values such as "High. sign." or " < "
#' should be coerced into numeric (default = `TRUE`).
#' 
#' @returns A list with 2 components:
#'     \itemize{
#'       \item \code{DataByPop}: a data.frame with 7 columns containing the full HWE output:
#'         \itemize{
#'           \item \code{Pop}: silly code
#'           \item \code{Locus}: locus name
#'           \item \code{PValue}: HWE p-value
#'           \item \code{SE}: standard error of the estimated p-value
#'           \item \code{WC Fis}: Weir and Cockerham Fis estimate
#'           \item \code{RH Fis}: Robertson and Hill Fis estimate
#'           \item \code{Steps}: number of genotypic matrices considered (exact) or number of switches (MCMC)
#'           }
#'       \item \code{SummaryPValues}: a matrix containing p-values for each locus (row) and pop (column) including the 
#'       overall p-value across pops and loci (Fisher's method)
#'       }
#'
#' @details
#' In `genepop` version 4.6 or below, `summaryAsNumeric = TRUE`, forces character output such as "High. sign." and to "0.000".
#' In `genepop` version 4.7 or avove, `summaryAsNumeric = TRUE`, forces character output such as " < 0.1555" and to "0.1555".
#' If MCMC, replace hard 0 p-values with `1 / (batches * iterations)` and calculate overall pops and overall loci via Fisher's method (ChiSqaure).
#' if Exact test, p-values overall pops or overall loci are pulled directly from `genepop` ("*.P") file
#'
#' @seealso [genepop::genepop-package()]
#'
#' @examples
#' \dontrun{
#' genepop_hwe <- read_genepop_hwe(file = "~/R/test.txt.P", sillyvec = sillyvec)
#' }
#' 
#' @export
read_genepop_hwe <- function(file, sillyvec = NULL, summaryAsNumeric = TRUE){

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
    HWdata=rbind(HWdata,cbind(Pop=rep(pops[i],ndiploci),reshape2::colsplit(
      sapply((popstart[i]+6):(popstart[i]+6+ndiploci-1),function(row){gsub(pattern="[[:blank:]]+",x=hwp[row],replacement="/",fixed=F)}),
      pattern="/",names=c("Locus","PValue","SE","WC Fis","RH Fis","Steps","Switches","Low"))[,1:6]))
  }
  
  # Break if Pop = 1
  
  if(npops == 1) {print(cat("Only 1 POP detected, so no summary SummaryPValues.\nAll p-values are pulled directly from the Genepop *.P file, none are derived in this function!\nHaploid loci are not included in the output data.frame."))
                  return(HWdata)
  }
  
  # Proceed as normal
  
  loci=gsub(reshape2::colsplit(hwp[grep("Locus ",hwp)],pattern="\"",names=c("Locus","Locus"))[,2],pattern="\"",replacement="")
  
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
    
    overallpops=reshape2::colsplit(hwp[grep("Locus ",hwp)+5+npops+4],pattern=" Prob :    ",names=c("Trash","P-value"))[,2]
    
    overallloci=reshape2::colsplit(hwp[grep("Pop : ",hwp)+5+ndiploci+4],pattern=" Prob :    ",names=c("Trash","P-value"))[,2]
    
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