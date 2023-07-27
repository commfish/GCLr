remove_ind_miss_loci <- function(sillyvec, loci = LocusControl$locusnames, proportion = 0.8, LocusCtl = LocusControl){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function removes individuals from "*.gcl" objects that have fewer non-missing loci than that specified by "proportion".
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   proportion - the cut-off proportion of the number of non-missing loci.
  #
  #   loci - optional vector of locus names to be considered when removing individuals for missing data; default is LocusControl$locusnames
  #          This argument is useful if you have loci in your markersuite that do not perform well in the lab and they will be dropped from the final dataset. 
  #   
  #   LocusCtl - an object created by [GCLr::create_locuscontrol()], (default = LocusControl) 
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Returns a tibble with 3 variables: 
  #                  SILLY_CODE <chr> = the silly with IDs removed 
  #                  IDs <list> = the IDs removed
  #                  is_empty <lgl> = were all IDs removed?
  #
  #    Assigns the ".gcl" objects back to workspace after removing individuals with missing loci.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  create_locuscontrol(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
  #  sillyvec <- c("SMCDO03", "SNEVA13")
  #  loki2r(sillyvec = sillyvec, username = "awbarclay", password = password)
  #
  #  naloci <- dimnames(SMCDO03.gcl)[[2]][100:199]
  #  empty.gcl <- SMCDO03.gcl
  #  empty.gcl[,naloci] <- NA
  #
  #  remove_ind_miss_loci(sillyvec = c(sillyvec, "empty"), loci = LocusControl$locusnames, proportion = 0.8)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!all(loci %in% LocusCtl$locusnames)){
    
    stop(paste0("The following `loci` were not found in `LocusControl`:\n", paste(setdiff(loci, LocusCtl$locusnames), collapse = "\n")))
    
  }
  
  nloci = length(loci)
  
  output <- lapply(sillyvec, function(silly){
    
    my.gcl <- get(paste0(silly, ".gcl"), pos = 1)
    
    tmp <- my.gcl %>% 
      dplyr::select(tidyselect::all_of(loci))  # subset for allele 1 of loci, because allele 2 defaults to NA for haploid
    
    IDsToRemove <- my.gcl %>% 
      dplyr::mutate(n_missing = rowSums(is.na(tmp))) %>% 
      dplyr::mutate(prop_loci = 1 - (n_missing/nloci)) %>% 
      dplyr::filter(prop_loci <= proportion) %>% 
      dplyr::pull(FK_FISH_ID)
    
    if(!purrr::is_empty(IDsToRemove)) {
      
      remove_ids(silly = silly, IDs = IDsToRemove)
      
    }

  }) %>% 
    dplyr::bind_rows()
  
  if(dim(output)[1] == 0){
    
    message("No individuals were removed")
    
  } else {
    
    message(paste0("A total of ", output$IDs_Removed %>% unlist() %>% length()), " individuals were removed from sillys in sillyvec.")
    
  }
  
  return(output)
  
}