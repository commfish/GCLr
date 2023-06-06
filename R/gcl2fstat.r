gcl2fstat <- function(sillyvec, loci, path, ncores = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   Write out an FSTAT input file from "*.gcl" objects.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - a character vector of locus names
  #   
  #   path - full file path to write out the FSTAT file with "\\" or "/" separator between folders 
  #           example: "V:\\Analysis\\2_Central\\Chinook\\Cook Inlet\\2019\\2019_UCI_Chinook_baseline_hap_data\\FSTAT\\FSTATfile.dat"
  #                   or "V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/FSTAT/FSTATfile.dat"
  #
  #   ncores - a numeric vector of length one indicating the number of cores to use
  #
  #   usat - logical vector of length 1, set to TRUE if loci are microsatellites
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   Writes out an FSTAT file to path.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #
  #   gcl2fstat(sillyvec = sillyvec31, loci = loci82, dir = "FSTAT", ncores = 8)
  # 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  start_time <- Sys.time()
  
  if(!exists("LocusControl")){
    
    stop("'LocusControl' not yet built.")
    
  }
  
  
  if(sum(is.na(match(loci, LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci,LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  nsillys <- length(sillyvec)
  
  nloci <- length(loci)
  
  ploidy <- LocusControl$ploidy[loci]
  
  alleles <- LocusControl$alleles[loci] %>% 
    dplyr::bind_rows(.id = "locus")
      
  maxchar <- max(nchar(alleles$allele))+1
    
  nalleles <- LocusControl %>% 
    dplyr::filter(locusnames%in%loci) %>% 
    dplyr::pull(nalleles) %>% 
    purrr::set_names(loci)
  
  my.gcl <- lapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)}) %>% 
    purrr::set_names(sillyvec)
  
  # Multicore loop
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  `%dopar%` <- foreach::`%dopar%`
  
  scores_all <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse")) %dopar% {
    
    new.gcl <- my.gcl[[silly]]
    
    pop_no <- match(silly, sillyvec)
    
    scores_names <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
      as.vector() 
    
    scores <- new.gcl %>% 
      dplyr::select(tidyselect::all_of(scores_names)) %>% 
      as.data.frame(stringsAsFactors = FALSE)
    
    lapply(loci, function(loc){
      
      variables <- c(loc, paste(loc, 1, sep = "."))
      
      my.alleles <- alleles %>% 
        dplyr::filter(locus==loc)
      
      scores %>%
        dplyr::select(tidyselect::all_of(variables)) %>% 
        dplyr::mutate(across(dplyr::everything(), .fns = ~factor(., levels = my.alleles$call))) %>% 
        dplyr::mutate(across(dplyr::everything(), .fns = as.numeric)) %>% 
        dplyr::mutate(across(dplyr::everything(), .fns = as.character)) %>%
        tidyr::replace_na(replace = list("0", "0") %>%
                            set_names(variables)) %>%
        dplyr::mutate(across(dplyr::everything(), .fns = ~stringr::str_pad(., width = maxchar, pad = "0", side = "left"))) %>% 
        tidyr::unite(col = !!rlang::as_name(loc), tidyselect::all_of(variables), sep = "")
      
    }) %>% 
      dplyr::bind_cols() %>% 
      purrr::set_names(loci) %>% 
      dplyr::mutate(Pop = pop_no) %>% 
      tidyr::unite("comb_scores", c(Pop, tidyselect::all_of(loci)), sep = " ") 
  
  } # end multicore loop
  
  parallel::stopCluster(cl)
  
  fstat <- paste(nsillys, nloci, max(nalleles), maxchar, sep = " ")
  
  fstat <- rbind(fstat, cbind(loci))
  
  fstat <- rbind(fstat,  scores_all %>% unlist() %>% cbind())
  
  write.table(x = fstat, file = path, row.names = FALSE, col.names = FALSE, quote = FALSE) 
  
  message("\nfstat file created\n", sep = '')
  
  Sys.time()-start_time
  
}