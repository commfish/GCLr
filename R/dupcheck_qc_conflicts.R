#' Investigate QC Samples With Conflicts
#' 
#'  @description
#'  This function checks samples with conflict rates > `conflict_rate` and runs [GCLr::dupcheck_among_sillys()] to see if they match any samples in `project_sillys`.(see details) 
#'  
#'  @param conflicts a tibble containing 3 variables:
#'  \itemize{
#'    \item \code{SillySource} original silly code and fish ID for each individual with conflicts
#'    \item \code{n} the number of conflicts among loci
#'    \item \code{p} the proportion of conflicts among loci
#'  }
#'  
#'  @param conflict_rate the conflict rate for determining which samples to investigate. This is assigned at the beginning of the QC script.
#'  @param project_sillys a character vector of silly codes included in the project
#'  @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = `LocusControl`)
#'  @param ncores the number of cores to run in parallel (default = [parallel::detectCores()])
#'  
#'  @seealso [GCLr::create_locuscontrol()]
#'  @seealso [GCLr::dupcheck_among_sillys()]
#'  
#'  @details
#'  The conflict rates used by this function are derived from Loki QC conflict reports. 
#'  Samples with high conflict rates may indicate catostrophic lab errors. 
#'  Checking these samples for duplicates in the project silly can help the QC analyzer determine if an error occured and how to fix the issue. 
#'  
#'  @note this function is only used in the QC script, so it's not exported to the package namespace
#'  
#'  @keywords invisible
dupcheck_qc_conflicts <- function(conflicts, conflict_rate, project_sillys, LocusCtl = LocusControl, ncores = parallel::detectCores()){
  
 # Filter for conflicts > conflict_rate
 conflicts_investigate <- conflicts %>%
   dplyr::filter(p > conflict_rate)
  
 # run the duplicate check, if necessary
 if (nrow(conflicts_investigate) == 0) {
  
  message(paste0("No individuals have > ", conflict_rate * 100, "% loci with conflicts between project and qc."))
    
  dup_check_results <- tibble::tibble(x = "Not applicable")
  
  } else{
    
    message(
      
      paste0("The following individuals have > ", conflict_rate * 100, "% loci with conflicts between project and qc:\n"),
      
      paste(conflicts_investigate$SillySource, conflicts_investigate$n, "conflicts", collapse = "\n")
      
      )
    
    # Loop through individuals to see if missing loci are an issue
    conflict_indiv <- sapply(conflicts_investigate$SillySource, function(silly_ind) {
      
      silly <- stringr::str_split(string = silly_ind, pattern = "_", simplify = TRUE)[, 1]
      
      ind <- stringr::str_split(string = silly_ind, pattern = "_", simplify = TRUE)[, 2]
      
      # qc fish lost in QA?
      if(ind %in% miss_loci_qc[[paste0(silly, "qc")]]){
        
        message(paste0("\n", silly, "qc_", ind," does not have at least 80% loci genotyped, not running DupCheck for this individual."))
        
      }  # print if any qc fish were removed due to missing genotypes
      
      # Project fish lost in QA
      if(ind %in% MissLoci[[silly]]){
        
        message(paste0("\n", silly, "_", ind, " does not have at least 80% loci genotyped, not running DupCheck for this individual."))
        
      }  # print if any project fish were removed due to missing genotypes
      
      paste(silly, ind, sep = "_")[!(ind %in% miss_loci_qc[[paste0(silly, "qc")]] | ind %in% MissLoci[[silly]])]  # Confirm qc fish and Project fish were not removed
      
    })  # silly_ind
    
    # If no more, stop
    
    if(is.null(conflict_indiv) | length(conflict_indiv) == 0){
      
      message("\nNo remaining high conflict individuals.")
      
      dup_check_results <- tibble::tibble(x = "Not applicable")
      
    } else{
      
      conflicts_investigate <- conflicts_investigate %>%
        dplyr::filter(SillySource %in% conflict_indiv)
      
      message("\nRunning dupcheck_among_sillys on these high conflict individuals, as they have at least 80% loci genotyped for Project and qc extractions.")
      
      message(paste(conflicts_investigate$SillySource, conflicts_investigate$n, "conflicts", collapse = "\n"))
      
      conflict_silly <- (stringr::str_split(string = conflicts_investigate$SillySource, pattern = "_", simplify = TRUE)[, 1]) %>%
        unique()
      
      KeySillyIDs <- lapply(conflict_silly, function(silly){
        
        tmp <- grep(pattern = silly, x = conflict_indiv, value = TRUE)
        
        sapply(tmp, function(ind){
          
          stringr::str_split(string = ind, pattern = "_", simplify = TRUE)[, 2]
          
        }, USE.NAMES = FALSE)
        
      }) %>% setNames(paste0(conflict_silly, "qc"))
      
      dup_check_results <- GCLr::dupcheck_among_sillys(KeySillys = paste0(conflict_silly, "qc"),
                                                       KeySillyIDs = KeySillyIDs,
                                                       BetweenSillys = project_sillys,
                                                       loci = loci,
                                                       minnonmissing = 0.6, 
                                                       minproportion = 0.9, 
                                                       ncores = ncores, 
                                                       plot.results = FALSE, 
                                                       LocusCtl = LocusCtl)
      
      if(is.null(dup_check_results)){
        
        dup_check_results <- tibble::tibble(x = "Not applicable")
        
      }
      
    }  # conflict_ind, post missing individuals
    
 }  # else, conflicts_to_investigate
  
  return(dup_check_results)
  
}