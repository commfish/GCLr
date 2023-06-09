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
#' 
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' loci <- GCLr::ex_baseline[,-c(1:5)] %>%
#'   names() %>%
#'   gsub(pattern = "*\\.1", x = ., replacement = "") %>%
#'   unique()
#'   
#' GCLr::pairwise_fst(sillyvec = sillyvec, loci = loci, inputfile = system.file("genepop", "ex_genepop.txt", package = "GCLr"), ncores = 4)
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