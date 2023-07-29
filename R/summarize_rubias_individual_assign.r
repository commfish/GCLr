#' Summarize `rubias` individual assignments.
#' 
#' This function will summarize `rubias` individual assignments.
#' 
#' @param rubias_output Output list object from \code{\link{run_rubias_mix}} or \code{\link{infer_mixture}}.
#' @param mixnames A character vector of mixture names to include in the individual assignment summary.
#' @param path Character vector of where to find output from each mixture as a .csv (created by \code{\link{run_rubias_mix}}).
#' 
#' @return A tibble with the following variables:
#'   - mixture_collection: The name of the mixture.
#'   - indiv: The individual (aka SillySource).
#'   - And a variable for each group in group_names containing the posterior means of group membership for each individual.
#' 
#' @examples
#' \dontrun{
#' path <- "C:/Users/awbarclay/Documents/Analysis/Sockeye/UCI_sockeye_2020_postseason/rubias/output"
#' rubias_output <- readRDS(file = "C:/Users/awbarclay/Documents/Analysis/Sockeye/UCI_sockeye_2020_postseason/output/rubias_output.rds")
#' mixnames <- c("DriftDW_20", "DriftCorr_20", "UpperSub_20", "Kasilof600ft_20", "WestKalgin_20", "Eastern_20", "General_North_20", "General_South_20")
#' summarize_rubias_individual_assign(rubias_output = rubias_output, mixnames = mixnames)
#' summarize_rubias_individual_assign(rubias_output = NULL, mixnames = mixnames, path = path)
#' }
#' 
#' @export
#' 
summarize_rubias_individual_assign <- function(rubias_output = NULL, mixnames = NULL, path = "rubias/output"){
  
  # Summarize rubias output files
  if(is.null(rubias_output)){
    
    if(is.null(mixnames)){
      
      stop("If rubias_output is NULL, then a vector of mixnames must be supplied for reading in the individual probabilities for each mixture.")
    
      }else{
        
      file_check <- sapply(mixnames, function(mix){
        
        !file.exists(paste0(path, "/", mix, "_indiv_posteriors.csv"))
        
      })
      
      if(sum(file_check) > 0){
        
        missing <- names(file_check[file_check])
        
        stop("No individual posterior files exist in the supplied path for the following mixtures: ", paste0(missing, collapse = ", "))
        
      }
      
     indiv_post <- lapply(mixnames, function(mix){
        
        readr::read_csv(file = paste0(path, "/", mix, "_indiv_posteriors.csv"), col_types = cols(
          mixture_collection = readr::col_character(),
          indiv = readr::col_character(),
          repunit = readr::col_character(),
          collection = readr::col_character(),
          PofZ = readr::col_double(),
          log_likelihood = readr::col_double(),
          z_score = readr::col_double(),
          n_non_miss_loci = readr::col_double(),
          n_miss_loci = readr::col_double()
        ))
        
      }) %>% dplyr::bind_rows()
      
      results <- indiv_post%>% 
        dplyr::mutate(repunit = factor(repunit, levels = unique(repunit))) %>% 
        dplyr::group_by(mixture_collection, indiv, repunit) %>% 
        dplyr::summarize(prob = sum(PofZ), .groups = "drop") %>% 
        tidyr::pivot_wider(names_from = repunit, values_from = prob)
      
    } 
    
    # Summarize rubias_output object
  } else{
    
    indiv_post <- rubias_output$indiv_posteriors
    
    mix_check <- !mixnames == indiv_post$mixture_collection %>% 
      unique() %>% 
      purrr::set_names(mixnames)
    
   if(sum(mix_check) > 0){
     
     stop("The folling mixnames are not included in the supplied rubias output object: ", paste0(names(mix_check[mix_check]), collapse = ", "))
     
   }
    
    results <- rubias_output$indiv_posteriors %>% 
      dplyr::filter(mixture_collection %in% mixnames) %>% 
      dplyr::mutate(repunit = factor(repunit, levels = unique(repunit))) %>% 
      dplyr::group_by(mixture_collection, indiv, repunit) %>% 
      dplyr::summarize(prob = sum(PofZ), .groups = "drop") %>% 
      tidyr::pivot_wider(names_from = repunit, values_from = prob)
    
  }
  
  return(results)
  
}
  
 