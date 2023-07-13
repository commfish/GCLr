#' @title Removes wrong species
#'
#' @description This function removes fish that have been identified, via genetic markers, as the wrong species. It relies on output from [find_alt_species()].
#'
#' @param AlternateSpeciesReport The object created by [find_alt_species()]
#' @param AlternateCutOff The minimum proportion of alternate alleles an individual must have in order to be considered the 'wrong species'. Default is 0.5
#' @param FailedCutOff The minimum proportion of failed loci an individual must have in order to be considered the 'wrong species'. Default is 0.5
#' @param NonmissingCutOff The minimum number of non-missing alternate loci an individual must have in order to be considered the 'wrong species'. Default is 3.
#'
#' @returns Returns a tibble showing which fish were removed and their stats (alternate and failed markers)
#'
#' @examples 
#' \dontrun{
#'  source(file = "Examples/qcExample.R")
#   wrong_spp <- GCLr::find_alt_species(sillyvec = sillys, species = "sockeye")
#   rm_spp <- GCLr::remove_alt_species(AlternateSpeciesReport = wrong_spp, AlternateCutOff = 0.5, FailedCutOff = 0.5, NonmissingCutOff = 3)
#' }
#' 
#' @export
remove_alt_species <- function(AlternateSpeciesReport, AlternateCutOff = 0.5, FailedCutOff = 0.5, NonmissingCutOff) {
  
  RemovedSpp <- AlternateSpeciesReport %>% # This is output from find_alt_species
    tidyr::separate(col = "silly_fish", into = c("silly", "ID"), sep = "\\_(?=[^\\_]+$)") %>% # Split out silly and fish
    dplyr::mutate(
      fate =
        dplyr::case_when(
          alternate >= AlternateCutOff & failure >= FailedCutOff & non_missing_alt >= NonmissingCutOff ~ "remove",
          TRUE ~ "keep"
        ) # If alternates AND failures AND # of non_missing_alts are greater than cutoffs, mark as remove
    ) %>%
    dplyr::filter(fate == "remove") # Just get fish marked to remove
  
  output <- lapply(RemovedSpp$silly %>% unique(), function(silly){
    
    IDsToRemove <-  RemovedSpp %>% 
      dplyr::filter(silly == !!silly) %>% 
      dplyr::pull(ID)
    
    if (!purrr::is_empty(IDsToRemove)) {
      
      GCLr::remove_ids(silly = silly, IDs = IDsToRemove)
      
    }else(data.frame())

  }) %>%
    dplyr::bind_rows()
  
  if (is_empty(output)) {
    
    warning("No fish were identified as the wrong species.")
    
  }
  
  return(output)
  
}

