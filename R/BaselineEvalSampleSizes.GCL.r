#' This function creates a tibble of sample sizes for creating baseline and mixture files for baseline evaluation tests.
#' The output from this function can be used by CreateRubiasBaselineEval.GCL() to write out mixture and baseline files for baseline evaluation tests in rubias.
#'
#' @param sillyvec A vector of silly codes without the ".gcl" extension.
#' @param group_names A character vector of group names with the length of max(groupvec).
#' @param groupvec A numeric vector indicating the group affiliation of each population in sillyvec.
#' @param mixsize The sample size of each test mixture.
#' @param scenarios A numeric vector of proportions to test for each group.
#' @param maxprop A numeric vector of length 1. This sets the maximum proportion of baseline individuals that can be 
#'                selected from each reporting group. For example, if maxprop = 0.5, the output tibble will only contain scenarios
#'                where sample sizes for the test_group do not exceed 50% of the fish in the baseline for that group.
#' @param seed An integer to set the seed for rmultinom(), so sample sizes are reproducible.
#'
#' @return The output of this function is a tibble with 4 variables: test_group, scenario, repunit, and samps.
#'         For each test_group and scenario, the number of rows will be length(group_names).
#'
#' @examples
#' BaselineEvalSampleSizes.GCL(sillyvec = c("KQUART06", "KQUART08", "KQUART10"), group_names = c("Upper_Susitna", "Chulitna", "Talkeetna", "Eastern_Susitna", "Deshka", "Yentna"), groupvec = c(1, 2, 3, 4, 5, 6), scenarios = round(seq(.01, 1, .01), 2), mixsize = 200, maxprop = 0.5, seed = 123)
  
BaselineEvalSampleSizes.GCL <- function(sillyvec, group_names, groupvec, mixsize, scenarios = round(seq(.01, 1, .01), 2), maxprop = 0.5, seed = 56){
  
  if(sum(str_detect(group_names, "\\W"))>0){
    
    stop("Special characters and spaces were detected in your group_names. 
          Using spaces and delimiters other than underscore in your group names may cause function errors later in your analysis.")
    
  }
    
  #Determine which scenarios to test for each group without removing maxprop of fish from the baseline.
  maxp <- silly_n.GCL(sillyvec) %>% 
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
      
      rs <- mixsize*p
      
      nr <- as.numeric(mixsize - rs)
      
      samps <- c(rs, rmultinom(n = 1, size = nr, prob = rep(1/(ngroups - 1), ngroups - 1))[ , 1]) %>% 
        purrr::set_names(c(g, setdiff(group_names, g)))
      
      tibble::tibble(test_group = g, scenario = p, repunit = names(samps), samps = samps)
      
    }) %>% dplyr::bind_rows()
    
  }) %>% dplyr::bind_rows()
    
}