locus_stats <- function(sillyvec, loci, ncores = 4, ...){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function creates a table of statistics containing obseved hetrozygosity (Ho), Fis, and Fst for each locus in loci.
  #   Note: as of 11/18/22 this function no longer writes out an FSTAT file using gcl2fstat, but instead calls on create_hierfstat_data to create 
  #         a hierfstat data object.Use gcl2FSTAT if you need to produce an FSTAT .dat file. 
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - a character vector of locus names
  #   
  #   ncores - a numeric vector of length one indicating the number of cores to use
  #
  #   ... - allows for previously used arguments 'fstatdir' and 'dir' to be supplied. 
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   A tibble is returned containing locusname, Ho, Fis, and Fst
  #
  #  if is.null(fstatdir) and exists(dir) "fstatfile.dat" is put into "dir"
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #
  #   locus_stats(sillyvec = sillyvec31, loci = loci82, ncores = 23)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  unused_args <- list(...)
  
  if (!length(unused_args) == 0){ 
    
    if(names(unused_args) %in% c("fstatdir", "dir")){
      
      warning("This function no longer writes out an FSTAT .dat file; therefore; the 'fstatdir' and 'dir' arguemnts are no longer used. 
            Use gcl2fstat to produce an FSTAT file if needed.")
      
      }
    
    }
  
  if(!exists("LocusControl")) {
    
    stop("'LocusControl' is required and not found, please create.")
    
  }
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  start.time <- Sys.time()
  
  ploidy <- LocusControl$ploidy[loci]
  
  dat <- GCLr::create_hierfstat_data(sillyvec = sillyvec, region = NULL, pop = seq_along(sillyvec), loci = loci, ncores = ncores)

  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores) #Start cluster  
  
  #Get variance components
  #Start parallel loop
  
  `%dopar%` <- foreach::`%dopar%`
  
  MyVC <- foreach::foreach(locus = loci, .packages = "hierfstat") %dopar% {
    
    diploid <- ploidy[locus]==2
    
    if(diploid){
      
      VC <- hierfstat::varcomp(dat[, c("pop", locus)], diploid = diploid)$overall
      
      table <- tibble::tibble(locus = !!locus, P = VC[1], I = VC[2], G = VC[3])
      
    }
    
    if(!diploid){
      
      VC <- hierfstat::varcomp(dat[, c("pop", locus)], diploid = diploid)$overall
      
      table <- tibble::tibble(locus = !!locus, P = VC[1], I = NA, G = VC[2])
      
    } 
    
    table
    
  } %>%  dplyr::bind_rows()
  
  parallel::stopCluster(cl)# Stop cluster
  
  #Summarize variance components
  MyTable0 <- MyVC %>% 
    dplyr::rowwise() %>%  
    dplyr::mutate(total = sum(c(P,I,G), na.rm = TRUE)) %>% 
    dplyr::mutate(P = P/total, I = I/total) %>% 
    dplyr::select(locus, P, I, total) %>% 
    dplyr::ungroup()
  
  Overall <- MyVC %>% 
    dplyr::summarize(P = sum(P, na.rm = TRUE), I = sum(I, na.rm = TRUE), total = sum(MyTable0$total, na.rm = TRUE)) %>% 
    dplyr::mutate(locus = "Overall", I = I/total, P = P/total)
  
  MyTable <- dplyr::bind_rows(MyTable0, Overall) 
  
  #Heterozygosities
  Hovec <- apply(hierfstat::basic.stats(dat[, c("pop", loci[ploidy==2])])$Ho, 1, mean, na.rm = TRUE)
  
  Ho <- tibble::tibble(locus = c(loci[ploidy==2], "Overall"), Ho = c(Hovec, mean(Hovec, na.rm = TRUE)))
  
  #Join varcomp summary with Ho  
  output <- MyTable %>% 
    dplyr::left_join(Ho, by = "locus") %>% 
    dplyr::mutate(Fis = I, Fst = P) %>% 
    dplyr::select(locus, Ho, Fis, Fst)
  
  print(Sys.time() - start.time)

  return(output)
  
}
