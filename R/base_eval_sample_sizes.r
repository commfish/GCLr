base_eval_sample_sizes <- function(sillyvec, group_names, groupvec, mixsize, scenarios = round(seq(.01, 1, .01), 2), maxprop = 0.5, seed = 56){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function creates a tibble of sample sizes for creating baseline and mixture files for baseline evaluation tests.
  #   
  #   The output from this function can be used by create_rubias_base_eval() to write out mixture and baseline files
  #   for baseline evaluation tests in rubias.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   group_names - character vector of group names the length of max(groupvec)
  #
  #   groupvec - numeric vector indicating the group affiliation of each pop in sillyvec
  #
  #   mixsize - the sample size of each test mixture
  #   
  #   scenarios - a numeric vector of proportions to test for each group. 
  # 
  #   maxprop - a numberic vector of length 1. The sets the maximum proportion of baseline individuals that can be 
  #             selected from each reporting group. For example, if maxprop = 0.5, the output tibble will only contain scenarios
  #             where sample sizes for the test_group do no exceed 50% of the fish in the baseline for that group.
  #
  #   seed - integer to set the seed for rmultinom(), so sample sizes are reproducible
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #  The output of this function is a tibble with 4 variables: test_group, scenario, repunit, and samps
  #  For each test_group and scenario, the number of rows will be length(group_names).
  #
  #  Here's an example output tibble for testing Upper_Susitna at 1% of the mixture..  
  #     test_group    scenario repunit         samps
  #     <chr>            <dbl> <chr>           <dbl>
  #     Upper_Susitna     0.01 Upper_Susitna       2  #Notice Upper_Susitna has 2 samples selected for a 200 sample mixuture (2/200 = 0.01)
  #     Upper_Susitna     0.01 Chulitna           39  #The remaining 198 samples are spread across the other reporting groups (repunits) 
  #     Upper_Susitna     0.01 Talkeetna          47
  #     Upper_Susitna     0.01 Eastern_Susitna    40
  #     Upper_Susitna     0.01 Deshka             35
  #     Upper_Susitna     0.01 Yentna             37
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  attach("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #  Final_Pops <- Final_Pops %>% mutate(group = factor(group, levels = unique(group)))
  #  base_eval_sample_sizes(sillyvec = Final_Pops$silly, group_names = Final_Pops$group %>% levels(), groupvec = Final_Pops$group %>% as.numeric(), scenarios = round(seq(.01, 1, .01), 2), mixsize = 200, maxprop = 0.5, seed = 123)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  if(sum(str_detect(group_names, "\\W"))>0){
    
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
      
      rs <- mixsize*p
      
      nr <- as.numeric(mixsize - rs)
      
      samps <- c(rs, rmultinom(n = 1, size = nr, prob = rep(1/(ngroups - 1), ngroups - 1))[ , 1]) %>% 
        purrr::set_names(c(g, setdiff(group_names, g)))
      
      tibble::tibble(test_group = g, scenario = p, repunit = names(samps), samps = samps)
      
    }) %>% dplyr::bind_rows()
    
  }) %>% dplyr::bind_rows()
    
}
  
  
 