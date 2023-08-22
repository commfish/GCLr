#' Sample Size by Locus
#'
#' This function takes a vector of sample names (sillyvec) and extracts data for specified loci from a 'LocusControl' object. It returns a data frame with the counts of non-missing values for each locus in each sample.
#'
#' @param sillyvec A vector of sample names to extract data from.
#' @param loci A character vector containing the names of loci for which data will be extracted. Defaults to the locus names in the 'LocusControl' object.
#' @param LocusCtl The 'LocusControl' object containing the data for all loci. Default is the global 'LocusControl' object.
#'
#' @return A data frame containing the counts of non-missing values for each locus in each sample. The output will have one row per sample and one column per locus, along with an additional 'silly' column indicating the sample name.
#'
#' @details
#' The function first checks if all specified loci in the 'loci' argument exist in the 'LocusControl' object. If any locus is not found, an error message is thrown.
#'
#' The function then iterates through each sample in 'sillyvec', retrieves the data for the specified loci from the 'LocusControl' object, and calculates the count of non-missing values for each locus. The results are combined into a single data frame.
#'
#' Note the output tibble from this function can be fed to Plot_SampleSizeByLocus.GCL to produce an interactive heatmap of the proportion of fish with scores for each locus and silly
#'
#' @seealso \code{\link{LocusControl}}
#'
#' @examples
#' \dontrun{
#' create_locuscontrol(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
#' sillyvec = c("SMCDO03", "SNEVA13")
#' password = "************"
#' loki2r(sillyvec = sillyvec, username = "awbarclay", password = password)
#' remove_ind_miss_loci(sillyvec = sillyvec)
#'   
#' sampn_by_locus(sillyvec = sillyvec, loci = LocusControl$locusnames)
#' }
#'
#' @export

sampn_by_locus <- function(sillyvec, loci = LocusControl$locusnames, LocusCtl = LocusControl){
  
  if (!all(loci %in% LocusCtl$locusnames)) {
    
    stop(paste0("'", setdiff(loci, LocusCtl$locusnames), "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
  }
  
  output <- lapply(sillyvec, function(silly){

    my.gcl <- get(paste0(silly, ".gcl"))
    
    my.gcl %>%
      dplyr::select(tidyselect::all_of(loci)) %>% 
      dplyr::summarize_all(function(x){sum(!is.na(x))}) %>% 
      dplyr::mutate(silly = silly) %>% 
      dplyr::select(silly, tidyselect::everything())
    
  }) %>% 
    dplyr::bind_rows()
  
  return(output)
  
}