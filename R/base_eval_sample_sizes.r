#' Get sample sizes for rubias baseline evaluation.
#' 
#' This function creates a tibble of sample sizes for creating baseline and mixture files for baseline evaluation tests.
#' 
#' @param sillyvec a character vector of silly codes without the ".gcl" extension (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
#' 
#' @param group_names character vector of group names the length of \code{max(groupvec)}
#'
#' @param groupvec a numeric vector indicating the group affiliation of each pop in `sillyvec`
#' 
#' @param mixsize a numeric vector of length 1; the sample size of each test mixture
#' 
#' @param scenarios a numeric vector of proportions to test for each group
#' 
#' @param  maxprop a numeric vector of length 1; the maximum proportion of baseline individuals to select from each reporting group (see details). 
#' 
#' @param seed integer to set the seed for [rmultinom()], so sample sizes are reproducible
#' 
#' @details The `maxprop` argument should be set to avoid oversampling populations, so the allele frequencies . For example, if `maxprop` = 0.5, the output tibble will only contain scenarios where sample sizes for the test_group do no exceed 50% of the fish in the baseline for that group.
#' 
#' @return a tibble with 4 variables: test_group, scenario, repunit, and samps. For each `test_group` and scenario, the number of rows will be \code{length(group_names)}.
#' 
#' @examples
#' 
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' group_names <- GCLr::ex_baseline$repunit %>% 
#'   unique()
#' 
#' groupvec <- GCLr::ex_baseline %>%
#'   dplyr::group_by(collection) %>%
#'   dplyr::filter(dplyr::row_number()==1) %>%
#'   dplyr::pull(repunit) %>%
#'   factor() %>%
#'   as.numeric()
#' 
#' GCLr::base_eval_sample_sizes(sillyvec = sillyvec, group_names = group_names, groupvec = groupvec, scenarios = round(seq(.01, 1, .01), 2), mixsize = 200, maxprop = 0.5, seed = 123)
#'
#' @export
base_eval_sample_sizes <- function(sillyvec, group_names, groupvec, mixsize, scenarios = round(seq(.01, 1, .01), 2), maxprop = 0.5, seed = 56){
  
  
  if(sum(stringr::str_detect(group_names, "\\W"))>0){
    
    stop("Special characters and spaces were detected in your group_names. 
          Using spaces and delimiters other than underscore in your group names may cause function errors later in your analysis.")
    
  }
    
  #Determine which scenarios to test for each group without removing maxprop of fish from the baseline.
  maxp <- silly_n(sillyvec) %>% 
    dplyr::mutate(groupvec = !!groupvec) %>% 
    dplyr::group_by(groupvec) %>% 
    dplyr::summarize(groupn = sum(as.numeric(n)), .groups = "drop_last") %>%
    dplyr::mutate(maxn = round(groupn*maxprop, 0), mixsize = mixsize) %>% 
    dplyr::mutate(maxprop = ifelse(round(maxn/!!mixsize, 1)>1, 1, round(maxn/mixsize, 1))) %>% 
    dplyr::pull(maxprop) %>% 
    purrr::set_names(group_names)
  
  ngroups <- length(group_names) 
  
  set.seed(seed)
  
  lapply(group_names, function(g){
    
    scenarios <- round(scenarios, digits = 2)
    
    gscn <- scenarios[scenarios <= maxp[g]]
    
    lapply(gscn, function(p){ 
      
      rs <- round(mixsize*p, digits = 0)
      
      nr <- as.numeric(mixsize - rs)
      
      samps <- c(rs, rmultinom(n = 1, size = nr, prob = rep(1/(ngroups - 1), ngroups - 1))[ , 1]) %>% 
        purrr::set_names(c(g, setdiff(group_names, g)))
      
      tibble::tibble(test_group = g, scenario = p, repunit = names(samps), samps = samps)
      
    }) %>% dplyr::bind_rows()
    
  }) %>% dplyr::bind_rows()
    
}