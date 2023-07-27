#' Remove Rows with NA Values for Each Element in a List
#'
#' This function removes any individuals from a SILLY.gcl object that were not genotyped for all markers in markersuite. It iterates through the list of individuals and removes the individuals with NA values in the specified columns for each element, storing the removed IDs in a tibble. The modified list is then returned.
#'
#' @param sillyvec a character vector of silly codes (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART09")) 
#' @param LocusCtl An optional LocusControl object (default: LocusControl) created by [GCLr::create_locuscontrol()] that specifies the columns to be checked for NA values. The function will use the \code{locusnames} attribute of the LocusCtl object to identify the columns.
#'
#' @return A tibble with two columns: \code{SILLY_CODE} contains the element names (SILLY_CODE) from \code{sillyvec}, and \code{IDs_Removed} contains the identifiers (IDs) of the rows removed for each element in the list.
#'
#' @details
#' The \code{remove_na_indv} function takes a list of elements, where each element is expected to be a character string representing a variable name. It then iterates through the list and processes each element separately. For each element, the function retrieves the corresponding data frame using the variable name and checks for NA values in the specified columns (taken from \code{LocusCtl$locusnames}). If a row contains any NA value in the specified columns, it is removed from the data frame. The function also creates a tibble that pairs the element name (SILLY_CODE) with the identifiers (IDs) of the removed rows (drop). Finally, all the tibbles are combined into one using \code{dplyr::bind_rows()} and returned as the output.
#'
#' @examples
#' \dontrun{
#' load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
#'  
#' removedInd <- remove_na_indv(sillyvec = sillyvec157)
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