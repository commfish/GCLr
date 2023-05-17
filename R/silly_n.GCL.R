#' Get sample size of .gcl objects
#' 
#' This function produces a tibble of .gcl object sample sizes
#' 
#' @param sillyvec A character vector of silly codes with associated *.gcl objects that you want sample sizes for.
#' 
#' @export
#' 
#' @return This function returns a two-column tibble of sillys (silly) and sample sizes (n)
#' 
#' @examples 
#' 
#' load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
#' 
#' silly_n.GCL(sillyvec = sillyvec157)
silly_n.GCL <- function(sillyvec){
  
  if(!require("pacman")) install.packages("pacman"); library(pacman); p_load(tidyverse) #Install packages, if not in library and then load them.
  
  lapply(sillyvec, function(silly){
    
    my.gcl <- get(paste0(silly, ".gcl"))
    
    if(!tibble::is_tibble(my.gcl)){
      
      stop("Make sure all of the sillys in sillyvec are in the new-style tibble format")
      
    }
    
    tibble::tibble(silly = silly, n = dim(my.gcl)[[1]])
    
  }) %>% 
    dplyr::bind_rows()
  
}