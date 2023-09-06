#' Remove Individuals with Missing Loci Data
#'
#' This function removes individuals from a data set based on the proportion of missing data in specified loci.
#'
#' @param sillyvec A character vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
#' @param loci A character vector specifying the names of loci (columns) in the datasets to be considered for missing data removal. Default is 'LocusControl$locusnames'.
#' @param proportion A numeric value (0 to 1) specifying the maximum proportion of missing data allowed for each individual in the specified loci. Default is 0.8.
#' @param LocusCtl An object of class 'LocusControl', containing locus names created by [GCLr::create_locuscontrol()]. Default is 'LocusControl'.
#'
#' @return A data frame containing the individuals removed from each dataset in 'sillyvec', along with their corresponding dataset names.
#'
#' @details
#' This function takes a vector of dataset names ('sillyvec') and removes individuals from each dataset based on the proportion of missing data in the specified loci. The 'loci' argument allows the user to specify the specific loci (columns) to be considered for missing data removal. By default, all loci in 'LocusControl' will be used.
#' 
#' The function calculates the proportion of missing data for each individual in the specified loci. If the proportion of missing data for an individual exceeds the 'proportion' argument, the individual is removed from the dataset. Only the individuals with a proportion of missing data less than or equal to 'proportion' will remain in the dataset.
#' 
#' Note that the function uses the 'get' function to access each dataset and relies on the 'dplyr' and 'purrr' packages for data manipulation and filtering.
#' 
#' If no individuals are removed from any of the datasets, the function will print a message indicating that no individuals were removed. Otherwise, it will print the total number of individuals removed from each dataset in 'sillyvec'.
#'
#' @examples
#' \dontrun{
#' create_locuscontrol(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
#' sillyvec <- c("SMCDO03", "SNEVA13")
#' loki2r(sillyvec = sillyvec, username = "awbarclay", password = password)
#' 
#' naloci <- dimnames(SMCDO03.gcl)[[2]][100:199]
#' empty.gcl <- SMCDO03.gcl
#' empty.gcl[,naloci] <- NA
#' 
#' remove_ind_miss_loci(sillyvec = c(sillyvec, "empty"), loci = LocusControl$locusnames, proportion = 0.8)
#' }
#'
#' @export

remove_ind_miss_loci <- function(sillyvec, loci = LocusControl$locusnames, proportion = 0.8, LocusCtl = LocusControl){
  
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