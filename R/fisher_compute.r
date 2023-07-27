#' @title Compute Fisher's Exact Test for Homogeneity of Allele Frequencies
#'
#' @description 
#' This function is a wrapper for [stats::fisher.test()] and runs a Fisher's Exact Test for homogeneity of allele frequencies.
#' The null hypothesis is that there are no significant differences in allele frequencies.
#'
#' @param freq A tibble of allele frequencies produced by [GCLr::calc_freq_pop()] for the collections (sillys) you want to test.
#' @param loci A character vector of locus names.
#' @param prec The precision of the output p-values (i.e., the number of significant digits; default = 4).
#'
#' @returns A tibble with 1 row with the following 3 columns:
#'     \itemize{
#'       \item \code{test_sillys}: a concatenated string of the collections (silly codes) being tested
#'       \item \code{overall}: overall p-value across all loci obtained using Fisher's method 
#'       \item \code{bylocus}: a nested list containing a tibble of locus-specific p-values
#'     }
#'
#' @details 
#' This function is called on by [GCLr::fishers_test()] to do multiple tests at once.
#' This function is used to test whether temporally or geographically similar baseline collections have significant differences in
#' allele frequencies across loci. If there are no significant differences between collections, then we typically pool collections into populations 
#' to improve estimates of allele frequencies.
#'
#' @examples
#' \dontrun{
#' load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
#' old2new_locuscontrol()
#' old2new_gcl(sillyvec = c("KKILL05","KKILL06"), save_old = FALSE)
#' freq <- calc_freq_pop(sillyvec = c("KKILL05","KKILL06"), loci = loci443)
#' fisher_compute(freq = freq, loci = loci443, prec = 4)
#' }
#'
#' @export
fisher_compute <- function(freq, loci, prec = 4) {
  
  # Getting p-values for each locus.
  pval <- sapply(loci, function(locus){
    
   my.freq <- freq %>% 
     dplyr::filter(locus == !!locus) %>% 
     dplyr::select(silly, allele_no, freq) 
   
   freq.mat <- matrix(my.freq$freq, ncol = max(my.freq$allele_no), dimnames = list(unique(my.freq$silly), unique(my.freq$allele_no)), byrow = TRUE)
   
   suppressWarnings(stats::fisher.test(freq.mat, workspace = 8000000, hybrid = TRUE))$p.value
    
  })
  
  pval[pval==0] <- .Machine$double.xmin #Replacing all hard zeros with the smallest nonnegative number possible in R.
  
  # Using Fisher's method to get overall p-value across loci
  overall <- stats::pchisq(q =- 2*sum(log(pval)), df = 2*length(loci), lower.tail = FALSE) %>% 
    round(prec)
  
  bylocus <- tibble::tibble(locus = names(pval), pval = pval %>% round(prec))
  
  test <- paste0(unique(freq$silly), collapse = ".")
  
  output <- tibble::tibble(test_sillys = test,  overall = !!overall, bylocus = list(!!bylocus))

  return(output)
  
}