pairwise_fst <- function(sillyvec, loci, inputfile, popnames = NULL, ncores = 4){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function generates a matrix of pairwise Fst values using the genepop package.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - a character vector of locus names
  # 
  #   inputfile - the file path of the genepop input file including .txt extension.
  #
  #   popnames - optional vector of population names corresponding to sillys in sillyvec to add as the dimnames of the output maxtrix. e.g., dimnames(ouput.matrix) <- list(popnames, popnames)
  #              If popnames is not supplied, sillyvec will be used as the dimnames e.g., dimnames(ouput.matrix) <- list(sillyvec, sillyvec)
  #   
  #   ncores - a numeric vector of length one indicating the number of cores to use
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  #  a matrix of pairwise fst values
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   source(paste0(path.expand("~/R/"), "Functions.R"))#GCL functions
  #
  #   load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #
  #   PWFST <- pairwise_fst(sillyvec = sillyvec31, loci = loci82, inputfile = "genepop/Susitna31pops82loci.txt", ncores = 8)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

  if(!file.exists(inputfile)){
    
    gcl2genepop(sillyvec = sillyvec, loci = loci, path = inputfile, ncores = ncores)
    
  }
  
  genepop::Fst(inputFile = inputfile, pairs = TRUE, outputFile = inputfile)
  
  if(is.null(popnames)){
    
    popnames = sillyvec
    
  }
  
 output <- read_genepop_pwfst(file = paste0(inputfile, ".MIG"), popnames = popnames)
 
 return(output)
  
}