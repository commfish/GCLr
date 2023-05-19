test_LD <- function(path, genepopFiles, dememorizations = 10000, batches = 100, iterations = 5000, ncores = parallel::detectCores() - 1){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function performs LD testing, using genepop::test_LD(), in parallel. 
  #   You can supply genepop files with n pops per file which are then tested for LD (these can be produced using gcl2genepop by supplying npops argument).  
  #   These DIS files are then combined into a single summary object containing p-valuse for each locus pair and population.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   path - the folder where the genepop files are located and where the results will be written
  #
  #   genepopFiles - the names of the genepop files located 
  #
  #   ncores - the number of cores for mulitcoring using doParallel and foreach. default = number of avaialble cores - 1
  #
  #   dememorizations - integer; length of dememorization step of Markov chain algorithm
  #
  #   batches	- integer; Number of batches
  #
  #   iterations - integer; Iterations per batch
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    When makeGenepops = TRUE, writes out a GENEPOP file with n pops per file.
  #    Writes out a GENEPOP DIS file with n pops per file. 
  #    Writes out a summarized GENEPOP DIS file with all results. 
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # source(paste0(path.expand("~/R/"), "Functions.R"))#GCL functions
  #
  # gcl2genepop(sillyvec, loci, path, VialNums = TRUE, usat = FALSE, ncores = 8, npops = 5) 
  #
  # my.files <- list.files(path = "GENEPOP", pattern = "(Pop)\\d+(to)\\d+(.gen.txt)") #The regular expression in the pattern argument finds the generic file names produced by gcl2genepop when npops is supplied.
  #
  # test_LD (genepopFiles = my.files, path = "GENEPOP", batches = 1, iterations = 1, ncores = 18)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  file.loci <- lapply(genepopFiles, function(file){ #Get locus names for output files and check to make sure all files contain the same loci
    
    my.file <- scan(paste0(path, "/", file), what = "")
    
    loci_start <- match(table = my.file, x = "format")+1
    
    loci_end <- match(table = my.file, x = "Pop")-1
    
    my.file[loci_start:loci_end]
    
  }) %>% unique()
  
  if(length(file.loci)==1){
    
    loci <- file.loci[[1]]
    
  }else{stop("Some genepop files contain different sets of loci.\nMake sure all files are created for the same loci.")}
  
  get.sillys <- function(file){ 
    
    my.file <- readr::read_lines(file = file)
    
    start <- match(table = my.file, x = "Pop")
    
    end <- length(my.file)
    
    tibble::as_tibble(my.file[start:end]) %>% 
      tidyr::separate(value, into = c("pop", NA), sep = " ,") %>% 
      tidyr::separate(pop, into = c("silly", NA), sep = "_") %>% 
      dplyr::filter(!silly == "Pop") %>% 
      dplyr::pull(silly) %>% 
      unique()
    
  }
  
  start_time <- Sys.time()
  
  original_directory <- getwd()
  
  dir.create(paste0(path, "/GenepopLD"))
   
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl) 
  
  LDout <- foreach::foreach(f = genepopFiles, .packages ="magrittr", .export = c("read_genepop_dis")) %dopar% { 
    
    setwd(original_directory) # set the original WD at beginning of each loop.
    
    new_directory <- paste0(path, "/GenepopLD/", f)
    
    dir.create(new_directory) # make a new directory for each genepop
    
    fs::file_move(path = paste0(path, "/", f), new_path = paste0(new_directory, "/", f)) # moving into newly created directory
    
    setwd(new_directory) # running from new directory to avoid overwriting files
    
    genepop::test_LD(inputFile = f, outputFile = paste0(f, ".DIS"), batches = batches, iterations = iterations) # run the test_LD function in parallel
    
    sillyvec = get.sillys(f) #This function extracts the silly codes from the genepop input file
    
    read_genepop_dis(file = paste0(f, ".DIS"), loci = loci, sillyvec = sillyvec)
    
    } %>% dplyr::bind_cols()#Combine results into a single DF
  
  parallel::stopCluster(cl) # kill cluster when complete
  
  dropcol <- c(dimnames(LDout)[[2]][grep(dimnames(LDout)[[2]], pattern = "(Locus)\\d\\.{3}\\d")][-c(1:2)],
               dimnames(LDout)[[2]][grep(dimnames(LDout)[[2]], pattern = "Overall")])#Duplicated columns to remove after combining outputs
  
  if(length(genepopFiles > 1)){ #The output will have different variables to drop depending on whether 1 or more than 1 genepop file is supplied.
    
    LDsummary <- LDout %>% 
      dplyr::mutate(Locus1 = `Locus1...1`, Locus2 = `Locus2...2`) %>% 
      dplyr::select(-dplyr::all_of(dropcol), -`Locus1...1`, -`Locus2...2`) %>% 
      dplyr::select(Locus1, Locus2, dplyr::everything())
    
  } else{ LDsummary <- LDout %>% 
    dplyr::select(-dplyr::all_of(dropcol)) %>% 
    dplyr::select(Locus1, Locus2, dplyr::everything())
  }
  
  # Combining locus pairs into a single column and put in long format.
  output <- LDsummary %>% 
    dplyr::mutate(Locus_pair = paste0(Locus1, "|", Locus2)) %>% 
    dplyr::select(Locus_pair, dplyr::everything()) %>% 
    dplyr::select(-Locus1, -Locus2) %>% 
    tidyr::pivot_longer(-Locus_pair, names_to = "Pop", values_to = "Pvalue") 
  
  Sys.time()-start_time
  
  return(output)
  
}  