create_bayes_base <- function(sillyvec, loci, dir, baseline_name, ncores = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates a baseline (.bse) file needed for `BAYES'.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - character vector of populations in the baseline
  #   loci - character vector of the loci you wish to include
  #   dir - character vector of where to save the ".bse" file
  #   baseline_name - character vector of what to name the baseline.bse, not including the .bse extension
  #   ncores - the number of cores to use in a foreach %dopar% loop. If the number of core exceeds the number on your device, then ncores defaults to detectCores() 
  #            
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns the fortran format of the baseline file - this object is needed for create_bayes_ctl
  #   Saves the baseline as a .bse file
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  # sillyvec <- Final_Pops$silly
  # loci <- loci82
  # base_fortran <- create_bayes_base(sillyvec = sillyvec, loci = loci, dir = getwd(), baseline_name = "SusitnaChinook31pops82loci", ncores = 8)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  if(!all(loci %in% LocusControl$locusnames)){
    
    stop(paste0("\n'", setdiff(loci, LocusControl$locusnames), "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
  }
  
  filename <- paste(dir, "\\", baseline_name, ".bse", sep = "")
  
  y <- GCLr::calc_freq_pop(sillyvec = sillyvec, loci = loci, ncores = ncores)
  
  bse0 <- y %>%
    dplyr::group_by(silly, locus) %>% 
    dplyr::mutate(n = sum(freq), silly = factor(silly, levels = sillyvec), locus = factor(locus, levels = loci)) %>% 
    dplyr::select(silly, locus, allele_no, n, freq) %>% 
    tidyr::pivot_wider(names_from = allele_no, values_from = freq, values_fill = 0) %>% 
    dplyr::mutate(silly = as.numeric(silly), locus = as.numeric(locus)) %>% 
    dplyr::arrange(silly, locus) %>% 
    dplyr::ungroup() 
  
  max_allele <- LocusControl$nalleles[loci] %>% 
    max()
  
  max_char <- bse0$n %>% 
    nchar() %>% 
    max()
  
  bse <- bse0 %>% 
    dplyr::mutate(dplyr::across(dplyr::everything(), ~stringr::str_pad(.x, width = 1 + max_char, side = "left"))) %>% 
    tidyr::unite(col = "file", dplyr::everything(), sep = "") %>% 
    dplyr::pull(file)
  
  base_fortran <- paste0("(",3 + max_allele, "(1X,I", max_char, "))")
  
  write.table(bse, filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(base_fortran)
  
}