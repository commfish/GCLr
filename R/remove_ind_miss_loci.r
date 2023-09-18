#' Remove Individuals with Missing Loci Data
#'
#' This function removes individuals from a data set based on the proportion of missing data in specified loci.
#'
#' @param sillyvec A character vector of silly codes without the ".gcl" extension (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
#' @param loci A character vector specifying the names of loci (columns) in the data sets to be considered for missing data removal. Default is 'LocusControl$locusnames'.
#' @param proportion A numeric value (0 to 1) specifying the maximum proportion of missing data allowed for each individual in the specified loci. Default is 0.8.
#' @param LocusCtl An object of class 'LocusControl', containing locus names created by [GCLr::create_locuscontrol()]. Default is 'LocusControl'.
#'
#' @return A data frame containing the individuals removed from each data set in 'sillyvec', along with their corresponding data set names.
#'
#' @details
#' This function takes a vector of data set (silly) names (`sillyvec`) and removes individuals from each silly based on the proportion of missing data in the specified loci. The `loci` argument allows the user to specify the specific loci (columns) to be considered for missing data removal. By default, all loci in `LocusControl` will be used.
#' 
#' The function calculates the proportion of missing data for each individual in the specified loci. If the proportion of missing data for an individual exceeds the `proportion` argument, the individual is removed from the silly. Only the individuals with a proportion of missing data less than or equal to `proportion` will remain in the silly.
#' 
#' Note that the function uses [base::get()] to access each data set and relies on the `dplyr` and `purrr` packages for data manipulation and filtering.
#' 
#' If no individuals are removed from any of the sillys, the function will print a message indicating that no individuals were removed. Otherwise, it will print the total number of individuals removed from each silly in `sillyvec`.
#'
#' @seealso [GCLr::remove_ids()]
#'
#' @examples
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' loci <- GCLr::ex_LocusControl$locusnames[-c(10, 12, 13, 32, 33, 97, 98)]
#' 
#' naloci <- dimnames(SURGOATM09.gcl)[[2]][100:201]
#' SURGOATM09.gcl[c(1, 5, 30), naloci] <- NA
#' assign("SURGOATM09.gcl", SURGOATM09.gcl, pos = -1)
#' 
#' GCLr::remove_ind_miss_loci(sillyvec = sillyvec, loci = loci, proportion = 0.8, LocusCtl = ex_LocusControl)
#' 
#' @export
remove_ind_miss_loci <- function(sillyvec, loci = LocusControl$locusnames, proportion = 0.8, LocusCtl = LocusControl){
  
  if(!all(loci %in% LocusCtl$locusnames)){
    
    stop(paste0("The following `loci` were not found in `LocusControl`:\n", paste(setdiff(loci, LocusCtl$locusnames), collapse = "\n")))
    
  }
  
  nloci <- length(loci)
  
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
      
      GCLr::remove_ids(silly = silly, IDs = IDsToRemove)
      
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