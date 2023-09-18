#' Get Sample Sizes from *.gcl Objects
#'
#' This function retrieves the sample sizes of *.gcl objects associated with given silly codes.
#'
#' @param sillyvec A character vector of silly codes with associated *.gcl objects for which you want to obtain sample sizes.
#'
#' @returns A two-column tibble containing the silly codes (silly) and their corresponding sample sizes (n).
#'
#' @examples
#'  
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' GCLr::silly_n(sillyvec = sillyvec)
#' 
#' @export
silly_n <- function(sillyvec) {
  
  lapply(sillyvec, function(silly){
    
    my.gcl <- get(paste0(silly, ".gcl"))
    
    if (!tibble::is_tibble(my.gcl)) {
      
      stop("Make sure all of the sillys in sillyvec are in the new-style tibble format")
      
    }
    
    tibble::tibble(silly = silly, n = dim(my.gcl)[[1]])
    
  }) %>% dplyr::bind_rows()
  
}