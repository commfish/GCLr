#' @title Perform Fisher's Exact Tests for Multiple Collections
#'
#' @description 
#' This function is a wrapper for [GCLr::fisher_compute()] to run Fisher's Exact Tests for homogeneity of allele frequencies
#' for multiple sets of collections. The null hypothesis is that there are no significant differences in allele frequencies.
#'
#' @param freq A tibble of allele frequencies produced by [GCLr::calc_freq_pop()] for the collections (sillys) you want to test.
#' @param loci A character vector of locus names.
#' @param tests A list of character vectors. Each element of the list represents a set of two or more collections (sillys) to be tested separately.
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()] (default = LocusControl).
#'
#' @returns A tibble with 1 row per `tests` with the following 3 columns:
#'     \itemize{
#'       \item \code{test_sillys}: a concatenated string of the collections (silly codes) being tested
#'       \item \code{overall}: overall p-value across all loci obtained using Fisher's method 
#'       \item \code{bylocus}: a nested list, named for each of `tests`, containing a tibble of locus-specific p-values:
#'         \itemize{
#'           \item \code{locus}: a character vector of locus names
#'           \item \code{pval}: a numeric vector of locus-specific p-values
#'         }
#'     }
#'
#' @details 
#' This function is used to test whether temporally or geographically similar baseline collections have significant differences in allele frequencies across loci.
#' If there are no significant differences between collections (overall p-value across loci > 0.01), then we typically pool collections into populations to improve 
#' estimates of allele frequencies.
#' 
#' @examples
#' \dontrun{
#' load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
#' old2new_locuscontrol()
#' old2new_gcl(sillyvec = c("KKILL05","KKILL06", "KFUNN05", "KFUNN06"), save_old = FALSE)
#' freq <- calc_freq_pop(sillyvec = c("KKILL05","KKILL06", "KFUNN05", "KFUNN06"), loci = loci443)
#' fishers_test(freq = freq, loci = loci443, tests = list(c("KKILL05","KKILL06"), c("KFUNN05", "KFUNN06")))
#' }
#'
#' @export
fishers_test <- function(freq, loci, tests, LocusCtl = LocusControl){

  if(!all(loci %in% LocusCtl$locusnames)){
    
    stop(paste0("'", setdiff(loci, LocusCtl$locusnames), "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
  }
  
  testnames <- sapply(tests, function(test){
    
    paste(test, collapse = ".")}
    
    )
  
  output <- lapply(tests, function(test){
    
    my.test.freq <- freq %>% 
      dplyr::filter(silly%in%test)
    
    GCLr::fisher_compute(freq = my.test.freq, loci = loci)
    
    }) %>% 
    dplyr::bind_rows() %>%
    dplyr::mutate(bylocus = bylocus %>% purrr::set_names(testnames))
    
  return(output)
  
}