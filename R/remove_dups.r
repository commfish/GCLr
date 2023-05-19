remove_dups <- function(dupcheck, remove_both = FALSE){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function removes duplicated "IDs" found with "dupcheck_within_silly"
  #  The default option is to remove the duplicated "IDs" with the highest number of missing loci for
  #  each duplicate pair (remove_both = FALSE). 
  #  If remove_both = FALSE and both IDs have the same number of missing loci, ID1 is removed.
  #  If remove_both = TRUE, both duplicated "IDs" are removed, use this option if you have paired data 
  #  that has become compromised by the duplicate pair (i.e. ASL, otolith data, etc.).
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   dupcheck - an object created by the function "dupcheck_within_silly". 
  #
  #   remove_both - TRUE/FALSE on whether to remove both duplicated "IDs" or just one of "IDs"
  #
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #    Returns a tibble with 3 variables: 
  #                  SILLY_CODE <chr> = the silly with IDs removed 
  #                  IDs <list> = the IDs removed
  #                  is_empty <lgl> = were all IDs removed?
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  password = "************"
  #  create_locuscontrol(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
  #  sillyvec = c("SMCDO03", "SNEVA13")
  #  loki2r(sillyvec = sillyvec, username = "awbarclay", password = password)
  #  remove_ind_miss_loci(sillyvec = sillyvec)
  #
  #  dupcheck <- dupcheck_within_silly(sillyvec = sillyvec, loci = LocusControl$locusnames, quantile = NULL, minproportion = 0.95, ncores = 8)
  #  removed_dups <- remove_dups(dupcheck)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(is_empty(dupcheck)){
    
    warning("Nothing removed. There are no duplicates to remove in dupcheck.", call. = FALSE)
    
    return(data.frame())
    
    }
  
  
  if(!tibble::is_tibble(dupcheck)){dupcheck <- dupcheck$report}
  
  dupcheck_names <-  c("silly", "ID1", "ID2", "Missing1", "Missing2", "proportion")
  
  if(!all(dupcheck_names %in% names(dupcheck))){
    
    stop(paste0("Nothing removed. Dupcheck must contain the folling variables: ", paste0(dupcheck_names, collapse = ", ")))
    
  }
  
  sillys <- dupcheck$silly %>% 
    unique()
  
  if(remove_both == TRUE) {
    
    to_remove <- dupcheck %>%
      tidyr::pivot_longer(
        cols = c("ID1", "ID2"),
        names_to = "ID",
        values_to = "remove"
      )
    
  } else {
    
    to_remove <- dupcheck %>%
      dplyr::mutate(
        remove = dplyr::case_when(
          Missing1 > Missing2 ~ ID1,
          Missing2 > Missing1 ~ ID2,
          Missing1 == Missing2 ~ ID1
        )
      )
    
  }
  
  output <- lapply(sillys, function(silly){
    
    remove <- to_remove %>% 
      dplyr::filter(silly == !!silly) %>% 
      dplyr::select(silly, removed_IDs = remove)
    
    remove_ids(silly = silly, IDs = remove$removed_IDs)
    
  }) %>% dplyr::bind_rows()
  
  return(output)
}