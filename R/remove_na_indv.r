#' Remove IDs with NA Values for all Loci
#'
#' This function removes any individuals from a SILLY.gcl object that were not genotyped for all Loci in markersuite.
#'
#' @param sillyvec a character vector of silly codes (e.g. `sillyvec <- c("KQUART06","KQUART08","KQUART09")`) 
#' @param LocusCtl An optional `LocusControl` object (default: `LocusControl`) created by [GCLr::create_locuscontrol()] that specifies the columns to be checked for NA values. The function will use the \code{locusnames} attribute of the LocusCtl object to identify the columns.
#'
#' @return A tibble with two columns: \code{SILLY_CODE} contains the element names (`SILLY_CODE`) from `sillyvec`, and `IDs_Removed` contains the identifiers (IDs) of the rows removed for each element in the list.
#'
#' @details
#' In Loki, NA genotypes mean that no genotyping was performed. 
#' This function checks each silly in `sillyvec` for individuals with NAs across all loci and removes them. 
#' This function should be run immediately after pulling genotypes using [GCLr::loki2r()], [GCLr::loki2r_gaps()], or [GCLr::loki2r_proj()] 
#'
#' @examples
#' \dontrun{
#' 
#' removedInd <- GCLr::remove_na_indv(sillyvec = sillyvec, LocusCtl = LocusControl)
#' 
#' }
#'
#' @export
remove_na_indv <- function(sillyvec, LocusCtl = LocusControl) {
  
  scores_cols <- LocusCtl$locusnames
  
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