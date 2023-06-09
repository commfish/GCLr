#' Read Pairwise Fst GENEPOP Output
#' 
#' This function reads in the pairwise FST half matrix from a GENEPOP "*.MIG" file and produces a full matrix with dimnames. The matrix from this function can be used to create a pairwise FST tree using [ape::nj()]
#' 
#' @param file the path to the GENEPOP .MIG file
#' 
#' @param popnames  a character vector of population names corresponding to the populations in the GENEPOP file used to create the "*.MIG" file
#' 
#' @param outfile the path including file name to save a tab delimited text file of the pairwise matrix
#' 
#' @return This function returns a named pairwise Fst matrix. If `!is.null(outfile)` the function produces a tab delimited text file of the pairwise matrix, otherwise no file is produced. 
#' 
#' @examples
#' GCLr::read_genepop_pwfst(file = system.file("genepop", "ex_genepop.txt.MIG", package = "GCLr"), popnames = paste0("Pop", seq(14)))
#' 
#' @export
read_genepop_pwfst <- function(file, popnames, outfile = NULL){

  npops <- length(popnames)
  
  f <- suppressWarnings(scan(file = file, skip = 3, nlines = ((npops-1)*npops/2), fill = TRUE, what = '', quiet = TRUE) %>% 
    as.numeric())
  
  ndist <- seq(((npops-1)*npops/2))
  
  f <- f[ndist]
  
  dims <- floor(sqrt(length(f) * 2))+1
  
  m <- matrix(NA, dims, dims)
  
  diag(m) <- 0
  
  m[upper.tri(m, diag = FALSE)] <- f
  
  m[lower.tri(m, diag = FALSE)] <- t(m)[lower.tri(t(m), diag = FALSE)]
  
  dimnames(m) <- list(popnames, popnames)
  
  if(!is.null(outfile)){
    
    write.table(m, file = outfile, sep="\t")
    
  } 
  
  m
  
}
