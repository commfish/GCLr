old2new_locuscontrol <- function (LocCtrl = LocusControl, save_old = FALSE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function converts the old style LocusControl object from a nested list to a tibble.
  # The old style LocusControl object can be saved as LocusControl_old before it is overwritten.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   LocCtrl - the locus control object
  #   overwrite - logical; whether you want the old LocusControl object overwritten without saving 
  #              (TRUE) or assign LocusControl to LocusControl_old before overwriting (FALSE) 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function assigns a converted LocusControl object to the current workspace.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  # 
  # old2new_locuscontrol(LocCtrl = LocusControl, save_old = TRUE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  if(tibble::is_tibble(LocCtrl)){
    
    stop("The LocusControl object is already a tibble.")
    
  }
  
  if(save_old){
    
    assign(paste0("LocusControl_old"), value = LocCtrl, pos = -1, envir = .GlobalEnv)
    
  }  # Saving old gcl
  
  # There was a bug with the old CombineLoci
  if(with(LocCtrl, length(locusnames) > length(Publishedlocusnames))){
    
    LocCtrl$Publishedlocusnames <- rep(NA, length(LocCtrl$locusnames))
    
    warning("Due to a bug with the old combine_loci() function, the length of Publishedlocusnames in the old LocusControl object is shorter than locusnames.
             To keep these variables the same length in the converted LocusControl tibble, Publishedlocusnames have been converted to NAs."
    )
  }
  
  locus_info <- with(LocCtrl, tibble::tibble(MarkerSuite, locusnames, Publishedlocusnames, nalleles, ploidy))
  
  alleles <- tibble::tibble(alleles = sapply(LocCtrl$locusnames, function(locus){
    
    a <- LocCtrl$alleles[[locus]]
    
    tibble::tibble(allele = seq(length(a)), call = a)
    
  }, simplify = FALSE)
  
  ) 
  
  LocCtrl_tidy = dplyr::bind_cols(locus_info, alleles)
  
  assign("LocusControl", value = LocCtrl_tidy, pos = -1, envir = .GlobalEnv)
  
  message("LocusControl converted from a list to a tibble")
  
}