#' Get Temporal Collections for Hierarchical ANOVA
#'
#' This function retrieves temporal collections from a vector of pooled sillys for hierarchical ANOVA. 
#' 
#' @note 
#' This function is intended for use while working within a baseline analysis workspace.  
#' Both pooled and unpooled sillys must be in your global environment or it will not work.
#' Single-silly pops and sillys with sample sizes < `min.samps` will not be included in the output.
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
#' @details
#' After pooling collections (sillys) from the same location and similar dates collected in different years (temporal collections) during baseline analysis, it's a good idea to see if allele frequencies are stable among sample years. 
#' But, first baseline analyzer must determine which temporal collections have sufficient sample size for accurate allele frequency estimates. 
#' With very large baselines, this can be a time consuming process. This function quickly gives you a table that can be used to create a `hierfstat` object in [GCLr::create_hierfstat_data()] to run a temporal Analysis of Variance in [GCLr::hierfstat_global_stats()] (see examples below). 
#' 
#' The optional region argument allows you to look at variance among regions, pops, and subpops.
#' 
#' @return a tibble with the following variables:
#'    \itemize{ 
#'      \item \code{silly}: non-pooled silly codes,
#'      \item \code{pop}: numbers indicating the population affiliation for each collection,
#'       \item \code{region}: numbers indicating the regional affiliation for each collection (Only included in output if `region` is supplied)
#'    }
#'
#' @examples
#' sillyvec_pool <- GCLr::base2gcl(GCLr::ex_baseline, unpool = FALSE)
#' 
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline, unpool = TRUE)
#' 
#' (temp_col <- GCLr::temporal_collections(sillyvec = sillyvec_pool, region = region, min.samps = 10, sep = ".")) #setting min.samps at 10 for this example, but really should be set at 50 or higher.
#' 
#' my.dat <- GCLr::create_hierfstat_data(sillyvec = temp_col$silly, pop = temp_col$pop, region = temp_col$region, loci = GCLr::ex_LocusControl$locusnames[-c(10, 12, 13, 32, 33)], LocusCtl = ex_LocusControl) %>%
#'   dplyr::mutate(region = factor(region), pop = factor(pop), spop = factor(spop))
#' 
#' GCLr::hierfstat_global_stats(levels = my.dat %>% dplyr::select(region, pop, spop), genotypes = my.dat %>% dplyr::select(-region, -pop, -spop), LocusCtl = ex_LocusControl)#Calculate variance components among regions, pops (pooled collections) within regions, sub pops (unpooled collections) within pops, and individuals within subpops.
#' 
#' GCLr::hierfstat_global_stats(levels = my.dat %>% dplyr::select(pop, spop), genotypes = my.dat %>% dplyr::select(-region, -pop, -spop), LocusCtl = ex_LocusControl)#Calculate variance components among pops (pooled collections), sub pops (unpooled collections) within pops, and individuals within subpops.
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