#' Create a hierfstat data object
#'
#' This function creates a \pkg{hierfstat} data object using the provided inputs. The \pkg{hierfstat} data object contains information about region (if supplied), population, sub-population numbers, and genotypes in single-column format. The object is used for computing hierarchical F statistics using functions such as [varcomp()] and [varcomp.glob()] from the [hierfstat] package.
#'
#' @param sillyvec a vector of silly codes without the ".gcl" extension
#' @param region optional; a numeric vector indicating the regional affiliation for each silly in `sillyvec`. Include this argument when computing hierarchical F statistics.
#' @param pop a numeric vector indicating the population affiliation for each silly in `sillyvec`. If there is only one silly per population in `sillyvec`, then this should be a sequence of numbers of length `sillyvec`.
#' @param loci a character vector of locus names.
#' @param ncores a numeric value indicating the number of cores to use.
#'
#' @return This function returns a [hierfstat] data object containing region (if supplied), population, sub-population numbers, and genotypes in single-column format.
#'
#' @examples
#' sillyvec <- c("KQUART06", "KQUART08", "KQUART10")
#' region <- c(1, 2, 3)
#' pop <- c(1, 1, 2)
#' loci <- c("locus1", "locus2", "locus3")
#' ncores <- 4
#' create_hierfstat_data(sillyvec, region, pop, loci, ncores)
#'
#' @import magrittr
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import foreach
#' @import doParallel
#' 
#' @aliases create_hierfstat_data.GCL
#'
#' @export

create_hierfstat_data <- function(sillyvec, region = NULL, pop,loci, ncores  = 4){

  if(!exists("LocusControl")){
    
    stop("'LocusControl' not yet built.")
    
  }
  
  if(sum(is.na(match(loci,LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci,LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  if(!all(sillyvec %in% stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))) {  # Do all sillys exist in the environment?
    
    missing_sillys <- setdiff(sillyvec, stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))
    
    stop(paste0("The following sillys are not in your environment:\n", paste0(missing_sillys, collapse = ", ")))
    
  }
  
  alleles <- LocusControl$alleles[loci] %>% 
    dplyr::bind_rows(.id = "locus")
  
  my.gcl <- lapply(sillyvec, function(silly){
    
    get(paste(silly,".gcl", sep = ""), pos = 1)
    
    }) %>% 
    purrr::set_names(sillyvec)
  
  # Multicore loop
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  output <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse")) %dopar% {
    
    new.gcl <- my.gcl[[silly]]
    
    spop_no <- match(silly, sillyvec)
    
    pop_no <- pop[spop_no]
    
    if(is.null(region)){
      
      region_no <- NA_integer_
      
    }else{region_no <- region[spop_no]}
    
    scores_names <- sapply(loci, function(locus) {
      
      c(locus, paste0(locus, ".1"))
      
      }) %>% as.vector() 
    
    scores <- new.gcl %>% 
      dplyr::select(tidyselect::all_of(scores_names)) %>% 
      as.data.frame(stringsAsFactors = FALSE)
    
    lapply(loci, function(loc){
      
      variables <- c(loc, paste(loc, 1, sep = "."))
      
      my.alleles <- alleles %>% 
        dplyr::filter(locus==loc)
      
      scores %>%
        dplyr::select(tidyselect::all_of(variables)) %>% 
        dplyr::mutate(across(dplyr::everything(), .fns = ~factor(., levels = my.alleles$call))) %>% 
        dplyr::mutate(across(dplyr::everything(), .fns = ~as.numeric(.))) %>%
        dplyr::mutate(across(dplyr::everything(), .fns = ~stringr::str_pad(., width = 2, side = "left", pad = "0"))) %>% 
        tidyr::unite(col = !!rlang::as_name(loc), tidyselect::all_of(variables), sep = "", na.rm = TRUE) %>% 
        dplyr::mutate(across(dplyr::everything(), .fns = ~as.numeric(.)))
      
      
    }) %>% 
      dplyr::bind_cols() %>% 
      purrr::set_names(loci) %>% 
      dplyr::mutate(region = region_no, pop = pop_no, spop = spop_no) 
  } %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(`region`, pop, spop, dplyr::everything())# end multicore loop
  
  parallel::stopCluster(cl)
  
  if(is.null(region)) {
    
    output <- output %>%
      dplyr::select(-`region`)
    
  }
  
  return(output)
  
}  

#' @rdname create_hierfstat_data
#' @export
create_hierfstat_data.GCL <- create_hierfstat_data

