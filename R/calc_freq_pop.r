calc_freq_pop <- function(sillyvec, loci = LocusControl$locusnames, ncores = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function gets the allele frequency for each locus in loci for each silly in sillyvec.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - vector of locus names; default is all LocusControl$locusnames.
  #
  #   ncores - the number of cores to use in a foreach %dopar% loop. If the nubmer of core exceeds the number on your device, then ncores defaults to detectCores()
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  A tibble with the following variables: silly, locus, allele_no (from LocusControl), freq (allele frequency), and proportion
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  #  old2new_locuscontrol()  
  #  old2new_gcl(sillyvec67)
  #
  #  Freq <- calc_freq_pop(sillyvec = sillyvec67, loci = loci413)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  start.time <- Sys.time() 
  
  if(!all(loci %in% LocusControl$locusnames)){
    
    stop(paste0("'", setdiff(loci, LocusControl$locusnames), "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
  }
  
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  all.gcl <- sapply(sillyvec, function(silly){get(paste0(silly, ".gcl"), pos = 1)}, simplify = FALSE)
  
  scores_cols <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
    as.vector() 
  
  alleles <- LocusControl$alleles[loci] %>% 
    dplyr::bind_rows(.id = "locus") %>% 
    dplyr::rename(allele_no = allele)
  
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  # Start parallel loop
  freqs <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse")) %dopar% {
    
    my.gcl <- all.gcl[[silly]]
    
    my.gcl %>% 
      dplyr::select(all_of(scores_cols)) %>% 
      tidyr::pivot_longer(cols = tidyselect::everything(), names_to = "locus", values_to = "allele") %>% 
      dplyr::mutate(locus = stringr::str_replace(string = locus, pattern = "\\.1$", replacement = "")) %>% 
      dplyr::count(locus, allele) %>% 
      dplyr::rename(freq = n) %>% 
      dplyr::right_join(alleles, by = c("locus" = "locus", "allele" = "call"), keep = TRUE) %>%
      tidyr::replace_na(list(freq = 0)) %>% 
      dplyr::mutate(silly = !!silly) %>% 
      dplyr::select(silly, locus = locus.y, allele_no, allele = call, freq) %>% 
      dplyr::arrange(locus, allele)
    
  } %>% 
    dplyr::bind_rows() 
  
  parallel::stopCluster(cl)  # End parallel loop
  
  output <- freqs %>% 
    dplyr::group_by(silly, locus) %>% 
    dplyr::mutate(total = sum(freq)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(proportion = freq/total) %>% 
    dplyr::select(-total)
  
  print(Sys.time() - start.time)
  
  return(output)
  
}