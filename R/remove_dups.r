#' Remove Duplicates based on Criteria
#'
#' This function removes duplicates from the input data frame based on specific criteria.
#'
#' @param dupcheck A data frame or tibble containing the information to check for duplicates. This object is produced by dupcheck_within_silly().
#' @param remove_both Logical. If set to TRUE, the function will remove duplicates from both "ID1" and "ID2" columns. If set to FALSE (default), the function will remove duplicates based on the "Missing1" and "Missing2" columns.
#'
#' @return A data frame containing the results after removing duplicates.
#' 
#' @details The `remove_dups` function takes an input data frame `dupcheck` and checks for duplicates based on the columns "silly", "ID1", "ID2", "Missing1", "Missing2", and "proportion". If there are no duplicates or if the input `dupcheck` is empty, the function will return an empty data frame with a warning message. The function will first check if `dupcheck` is a tibble; if not, it assumes the data is stored in a column named "report" and extracts the tibble from that column.
#' 
#' If the `remove_both` argument is set to TRUE, the function removes duplicates from both "ID1" and "ID2" columns using tidyr::pivot_longer. Otherwise, it calculates which duplicates to remove based on the "Missing1" and "Missing2" columns using dplyr::mutate and dplyr::case_when.
#' 
#' The function then processes the data for each unique value of "silly" and removes the corresponding duplicates using the `remove_ids` function.
#'
#' @examples
#' \dontrun{
#' password = "************"
#' create_locuscontrol(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
#' sillyvec = c("SMCDO03", "SNEVA13")
#' loki2r(sillyvec = sillyvec, username = "awbarclay", password = password)
#' remove_ind_miss_loci(sillyvec = sillyvec)
#' 
#' dupcheck <- dupcheck_within_silly(sillyvec = sillyvec, loci = LocusControl$locusnames, quantile = NULL, minproportion = 0.95, ncores = 8)
#' removed_dups <- remove_dups(dupcheck)
#' }
#' 
#' @seealso
#' \code{\link{remove_ids}}
#' 
#' @export
remove_dups <- function(dupcheck, remove_both = FALSE){
  
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