#' Create Pairwise Fst Matrix
#'
#' This function generates a matrix of pairwise Fst values using the \pkg{genepop} package.
#'
#' @param sillyvec a vector of silly codes without the ".gcl" extension (e.g., sillyvec <- c("KQUART06", "KQUART08", "KQUART10"))
#' 
#' @param loci a character vector of locus names
#' 
#' @param inputfile the file path of the \pkg{genepop} input file including the .txt extension
#' 
#' @param popnames optional vector of population names corresponding to `sillyvec` to add as the dimnames of the output matrix. If `popnames` is not supplied, `sillyvec` will be used as the dimnames.
#' 
#' @param ncores a numeric vector of length one indicating the number of cores to use
#'
#' @return a matrix of pairwise Fst values.
#'
#' @examples
#' source(paste0(path.expand("~/R/"), "Functions.R")) # GCL functions
#' load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
#' PWFST <- pairwise_fst(sillyvec = sillyvec31, loci = loci82, inputfile = "genepop/Susitna31pops82loci.txt", ncores = 8)
#'
#' @aliases PW_FST.GCL
#'
#' @export

pairwise_fst <- function(sillyvec, loci, inputfile, popnames = NULL, ncores = 4){
 
  if(!file.exists(inputfile)){
    
    GCLr::gcl2genepop(sillyvec = sillyvec, loci = loci, path = inputfile, ncores = ncores)
    
  }
  
  genepop::Fst(inputFile = inputfile, pairs = TRUE, outputFile = inputfile)
  
  if(is.null(popnames)){
    
    popnames = sillyvec
    
  }
  
 output <- GCLr::read_genepop_pwfst(file = paste0(inputfile, ".MIG"), popnames = popnames)
 
 return(output)
  
}
#' @rdname pairwise_fst
#' @export
PW_FST.GCL <- pairwise_fst  