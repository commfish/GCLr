#' Convert Old LocusControl to New LocusControl Tibble
#'
#' This function converts an old LocusControl object to a new LocusControl tibble
#' to accommodate changes in the data structure. The old LocusControl object is 
#' expected to be a list, and the new LocusControl tibble will have a more 
#' standardized and tidy format.
#'
#' @param LocCtrl An old LocusControl object (list) to be converted. If the 
#'   object is already a tibble, an error will be thrown.
#'
#' @param save_old Logical, indicating whether to save the old LocusControl 
#'   object as "LocusControl_old" in the global environment. Default is FALSE.
#'
#' @return The function does not return a value explicitly, but it creates a new 
#'   LocusControl tibble with updated format and saves it in the global 
#'   environment as "LocusControl".
#'
#' @details
#' This function takes an old LocusControl object and converts it to a new 
#' LocusControl tibble, which consists of two parts. The first part, locus_info, 
#' contains the following columns:
#'
#' \describe{
#'   \item{MarkerSuite}{The name of the MarkerSuite.}
#'   \item{locusnames}{A character vector of locus names.}
#'   \item{Publishedlocusnames}{A character vector of published locus names. 
#'     If the old LocusControl object has a shorter Publishedlocusnames vector 
#'     than locusnames, NA values will be used to make them the same length.}
#'   \item{nalleles}{The number of alleles for each locus.}
#'   \item{ploidy}{The ploidy level for each locus.}
#' }
#'
#' The second part, alleles, contains the alleles information for each locus. 
#' It has two columns:
#'
#' \describe{
#'   \item{alleles}{A list column containing the allele numbers for each locus.}
#'   \item{call}{A list column containing the alleles calls for each locus.}
#' }
#'
#' Note that this function will save the old LocusControl object as 
#' "LocusControl_old" in the global environment if the \code{save_old} argument 
#' is set to TRUE.
#'
#' @examples
#' \dontrun{
#' load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
#'
#' GCLr::old2new_locuscontrol(LocCtrl = LocusControl, save_old = TRUE)
#' }
#'
#' @seealso
#' \code{\link{LocusControl}}, \code{\link{combine_loci}}
#'
#' @export
old2new_locuscontrol <- function (LocCtrl = LocusControl, save_old = FALSE){
  
  if(tibble::is_tibble(LocCtrl)){
    
    stop("The LocusControl object is already a tibble.")
    
  }
  
  if(save_old){
    
    assign(paste0("LocusControl_old"), value = LocCtrl, pos = -1, envir = .GlobalEnv)
    
  }  # Saving old gcl
  
  # There was a bug with the old CombineLoci
  if(with(LocCtrl, length(locusnames) > length(Publishedlocusnames))){
    
    LocCtrl$Publishedlocusnames <- rep(NA, length(LocCtrl$locusnames))
    
    warning("Due to a bug with the old combine_loci() function, the length of Publishedlocusnames in the old LocusControl object is shorter than locusnames.
             To keep these variables the same length in the converted LocusControl tibble, Publishedlocusnames have been converted to NAs."
    )
  }
  
  locus_info <- with(LocCtrl, tibble::tibble(MarkerSuite, locusnames, Publishedlocusnames, nalleles, ploidy))
  
  alleles <- tibble::tibble(alleles = sapply(LocCtrl$locusnames, function(locus){
    
    a <- LocCtrl$alleles[[locus]]
    
    tibble::tibble(allele = seq(length(a)), call = a)
    
  }, simplify = FALSE)
  
  ) 
  
  LocCtrl_tidy = dplyr::bind_cols(locus_info, alleles)
  
  assign("LocusControl", value = LocCtrl_tidy, pos = -1, envir = .GlobalEnv)
  
  message("LocusControl converted from a list to a tibble")
  
}