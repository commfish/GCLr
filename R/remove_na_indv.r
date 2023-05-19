remove_na_indv <- function(sillyvec) {
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function removes any individuals from a SILLY.gcl object that were not genotyped
  #   for all markers in markersuite.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a character vector of silly codes (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART09")) 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    1) this function will modify silly.gcls in your environment to remove NA fish
  #    2) output is tibble of the fish IDs removed by silly
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  # 
  #   removedInd <- remove_na_indv(sillyvec = sillyvec157)
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function requires a locus control object to subset out the scores columns.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!exists("LocusControl")){
    
    stop("'LocusControl' not yet built.")
    
  }
  
  
  scores_cols <- LocusControl$locusnames
  
  na.individuals.removed <- lapply(sillyvec, function(silly) {
    
    my.gcl <- get(paste0(silly, ".gcl"))
    
    drop <- my.gcl$FK_FISH_ID[apply(my.gcl[ , scores_cols], 1, function(ind) {
      
      sum(is.na(ind)) > 0 
      
    })]
    
    assign(
      x = paste0(silly, ".gcl"),
      value = my.gcl %>%
        filter(!FK_FISH_ID %in% drop),
      pos = -1,
      envir = .GlobalEnv
    )
    
    tibble::tibble(SILLY_CODE = silly, IDs_Removed = drop)
    
  }) %>%
    dplyr::bind_rows()
  
  return(na.individuals.removed)
  
}