#' Get Temporal Collections for Hierarchical ANOVA
#'
#' This function retrieves temporal collections from a vector of pooled sillys for hierarchical ANOVA. 
#' 
#' @note Single-silly pops and sillys with sample sizes < `min.samps` will not be included in the output.
#'
#' @param sillyvec a vector of pooled silly codes (a.k.a. pooled pops) without the ".gcl" extension
#' 
#' @param region an optional tibble containing two numeric variables: 
#'             \itemize{
#'               \item \code{pop_no}: the position number for each population in sillyvec
#'               \item \code{region}: the region number for each population in sillyvec
#'               }
#'               
#' @param min.samps the minimum number of samples for including a collection
#' 
#' @param sep the separator used when pooling collections into populations; default is "."
#'
#' @return a tibble with the following variables:
#'    \itemize{ 
#'      \item \code{silly}: non-pooled silly codes,
#'      \item \code{spop}: numbers indicating the population affiliation for each collection,
#'       \item \code{sregion}: numbers indicating the regional affiliation for each collection
#'    }
#'
#' @examples
#' load("V:/Analysis/2_Central/Coho/Cook Inlet/2019/2019_Cook_Inlet_coho_baseline/2019_Cook_Inlet_coho_baseline_new.Rdata")
#' source(paste0(path.expand("~/R/"), "Functions.R"))#GCL functions
#' sillyvec <- Final_Pops$collection
#' region <- Final_Pops %>% select(pop = order, region = region)
#' temporal_collections(sillyvec = sillyvec, region = region, min.samps = 50, sep = ".")
#'
#' @aliases temporal_collections.GCL
#' 
#' @export 

temporal_collections <- function(sillyvec, region = NULL, min.samps = 50, sep = "."){
  
  if(!all(sillyvec %in% stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))) {  # Do all sillys exist in the environment?
    
    missing_sillys <- setdiff(sillyvec, stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))
    
    stop(paste0("The following sillys are not in your environment:\n", paste0(missing_sillys, collapse = ", ")))
    
  }
  
  tc <- lapply(1:length(sillyvec), function(pop){
    
    silly <- sillyvec[pop]
    
    tcol0 <- base::strsplit(silly, split = paste0("\\", sep))[[1]]
    
    tcol <- GCLr::silly_n(tcol0) %>% 
      dplyr::filter(n >= min.samps) %>% 
      dplyr::pull(silly)
    
    if(length(tcol) < 2){
      
      tibble::tibble(pop = pop, silly = NA_character_)
      
    }else(tibble::tibble(pop = pop, silly = tcol))
    
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(!is.na(silly)) %>% 
    dplyr::select(silly, pop)
  
  if(!is.null(region)){
    
    tc <- tc %>% 
      dplyr::left_join(region, by = c("pop"="pop_no")) %>% 
      dplyr::mutate(region = factor(region, levels = unique(region)) %>% as.numeric())
    
  }
  
  return(tc)
  
}

#' @rdname temporal_collections
#' @export
temporal_collections.GCL <- temporal_collections  