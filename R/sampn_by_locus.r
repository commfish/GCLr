sampn_by_locus <- function(sillyvec, loci = LocusControl$locusnames){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function creates a tibble of sample sizes by locus.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - vector of locus names.
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  #   A tibble of sample sizes by locus for each silly.
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  create_locuscontrol(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
  #  sillyvec = c("SMCDO03", "SNEVA13")
  #  password = "************"
  #  loki2r(sillyvec = sillyvec, username = "awbarclay", password = password)
  #  remove_ind_miss_loci(sillyvec = sillyvec)
  #
  #  sampn_by_locus(sillyvec = sillyvec, loci = LocusControl$locusnames)
  # 
  #  Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  The output tibble from this function can be fed to Plot_SampleSizeByLocus.GCL to produce an interactive heatmap 
  #  of the proportion of fish with scores for each locus and silly
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!all(loci %in% LocusControl$locusnames)){
    
    stop(paste0("'", setdiff(loci, LocusControl$locusnames), "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
  }
  
  
  nsilly <- length(sillyvec)
  
  output <- lapply(sillyvec, function(silly){
    
    my.gcl <- get(paste0(silly, ".gcl"))
    
    my.gcl %>%
      dplyr::select(tidyselect::all_of(loci)) %>% 
      dplyr::summarize_all(function(x){sum(!is.na(x))}) %>% 
      dplyr::mutate(silly = silly) %>% 
      dplyr::select(silly, tidyselect::everything())
    
  }) %>% 
    dplyr::bind_rows()
  
  return(output)
  
}