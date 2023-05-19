read_genepop_pwfst <- function(file, popnames, outfile = NULL){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function reads in the pairwise FST half matrix from a GENEPOP "*.MIG" file and produces a full matrix with dimnames.
  # The matrix from this function can be used to create a pairwise FST tree using ape::nj() 
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   file - the path to the GENEPOP .MIG 
  #   popnames - a character vector of population names corresponding to the populations in the GENEPOP file used to create the "*.MIG" file.
  #   outfile - the path including file name to save a tab delimited text file of the pairwise matrix. 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function returns a named pariwise matrix of FSTs.
  #   If !is.null(outfile) the function produces a tab delimited text file of the pairwise matrix, otherwise no file is produced. 
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # attach("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #
  # file <- "V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/GENEPOP/Sustina36pops82loci.MIG"
  # 
  # read_genepop_pwfst(file = file, popnames = Pooling_1_pops$location)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
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
