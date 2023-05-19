summarize_rubias_base_eval <- function(mixvec, sample_sizes, method = c("MCMC", "PB", "both")[1], group_names = NULL, group_names_new = NULL,
                                               groupvec = NULL, groupvec_new = NULL, path = "rubias/output", alpha = 0.1, 
                                               burn_in = 5000,  threshold = 5e-7, ncores = 4){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function is a wrapper for custom_comb_rubias_output(). The function creates a tibble of sample sizes for creating baseline and mixture files for baseline evaluation tests.
  #   
  #   The output from this function can be used by create_rubias_base_eval() to write out mixture and baseline files
  #   for baseline evaluation tests in rubias.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   mixvec - a vector of test mixture names
  # 
  #   sample_sizes - a tibble produced by base_eval_sample_sizes() containing the following variables: test_group, scenario, repunit, and samps
  #
  #   method - a character string indicating the rubias output to summarize.  Select one of three choices: "MCMC" (summarize MCMC output), "PB" (summarize bias corrected output), "both" (summarize both outputs)
  #
  #   see custom_comb_rubias_output() documentation for the remaining inputs.
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #  The output of this function is a list containing two elements: 1) estimates and 2) summary_stats.
  #  
  #  estimates - a tibble containing the following variables: test_group, scenario, repunit, mean, sd, lo5CI (5% CI), hi95CI (95% CI), P=0, true_proportion (actual proportion in mixture), total_samples (mixture size), and method (PB or MCMC)
  #  
  #  summary_stats - a tibble containing mixture summary statistics for each group tested: method (PB or MCMC), test_group, RMSE (root mean squared error), Mean_Bias (average bias among all test mixtures for a given test_group),
  #                  90%_within (the maximum distance from the true proportion where 90% of point estimates occurred), and Within_Interval (the proportion of tests where the credibility interval contained the true proportion).
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #  require(tidyverse)
  #  tests <- sample_sizes %>% group_by(test_group, scenario) %>% summarize(test_group = test_group %>% unique(), scenario = scenario %>% unique(), .groups = "drop_last")#Total of 510 tests  #
  #  mixvec <- tests %>% unite(col = "mixvec", test_group, scenario, sep ="_" ) %>% pull()
  #  path <-  "V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/rubias/output/3groups"
  # 
  #  summarize_rubias_base_eval (mixvec = mixvec, sample_sizes = sample_sizes, method = "both", group_names = NULL, group_names_new = NULL, groupvec = NULL, groupvec_new = NULL, path = path, alpha = 0.1,
  #                                  burn_in = 5000, threshold = 5e-7, ncores = 8)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  start_time <- Sys.time()
  

  # Get summary estimates
  ## MCMC
  if(method == "MCMC"){
    
    estimates <- custom_comb_rubias_output(rubias_output = NULL, mixvec = mixvec, group_names = group_names, group_names_new = group_names_new,
                                              groupvec = groupvec, groupvec_new = groupvec_new, path = path, alpha = alpha, 
                                              burn_in = burn_in, bias_corr = FALSE, threshold = threshold, plot_trace = FALSE, ncores = ncores) %>% 
      dplyr::mutate(method = method)
  
  }
  
  ## PB (bias correction)
  if(method == "PB"){
    
    
    estimates <- custom_comb_rubias_output(rubias_output = NULL, mixvec = mixvec, group_names = group_names, group_names_new = group_names_new,
                                              groupvec = groupvec, groupvec_new = groupvec_new, path = path, alpha = alpha, 
                                              burn_in = burn_in, bias_corr = TRUE, threshold = threshold, plot_trace = FALSE, ncores = ncores) %>% 
      dplyr::mutate(method = method)
    
  }
  
  ## MCMC & PB
  if(method == "both"){
    
    cl <- parallel::makePSOCKcluster(ncores)
    
    doParallel::registerDoParallel(cl, cores = ncores)  
    
    parallel::clusterExport(cl = cl, varlist = "custom_comb_rubias_output") # For some reason custom_comb_rubias_output() couldn't be found within the foreach loop, this solves that issue.
    
    estimates <- foreach::foreach(bias_corr = c(TRUE, FALSE), .packages = "tidyverse", .export = "custom_comb_rubias_output") %dopar% {
      
      custom_comb_rubias_output(rubias_output = NULL, mixvec = mixvec, group_names = group_names, group_names_new = group_names_new,
                                   groupvec = groupvec, groupvec_new = groupvec_new, path = path, alpha = alpha, 
                                   burn_in = burn_in, bias_corr = bias_corr, threshold = threshold, plot_trace = FALSE, ncores = ncores) %>% 
        dplyr::mutate(method = if(bias_corr){"PB"}else{"MCMC"})
      
    } %>% dplyr::bind_rows()
    
    parallel::stopCluster(cl)
    
  }
  
  # Create estimates output summary
  estimates_out <- estimates %>% 
    tidyr::separate(mixture_collection, into = c("test_group", "scenario"), sep = "\\_(?=[^\\_]+$)", remove = TRUE) %>% #The regular expression means the last instance of "_", just in case someone has a group name with an underscore.
    dplyr::mutate(scenario = as.numeric(scenario)) %>% 
    dplyr::left_join(sample_sizes, by = c("repunit" = "repunit", "scenario" = "scenario", "test_group" = "test_group")) %>% 
    dplyr::group_by(test_group, scenario, method) %>% 
    dplyr::mutate(total_samps = sum(samps)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(repunit = factor(repunit, levels = unique(repunit)), true_proportion = samps/total_samps, lo5CI = `5%`, hi95CI = `95%`) %>% 
    dplyr::select(test_group, scenario, repunit, mean, sd, lo5CI, hi95CI, `P=0`, true_proportion, total_samps, method)
  
  # Calculate the baseline eval summary statistics RMSE =  root mean squared error, 'Within_10%' =  proportion of estimates within 10% of true
  summary_stats <- estimates_out %>% 
    dplyr::filter(test_group==repunit) %>% 
    dplyr::group_by(method, test_group) %>%
    dplyr::mutate(Bias = mean-true_proportion, Abs_Bias = abs(mean-true_proportion), Bias_squared = (mean-true_proportion)^2, CI_width = hi95CI-lo5CI,
           lo5CI_within = true_proportion-lo5CI, hi95CI_within = hi95CI-true_proportion, Within_CI = dplyr::case_when(true_proportion >= lo5CI & true_proportion <= hi95CI ~ TRUE, TRUE~FALSE)) %>%
    dplyr::summarise(RMSE = sqrt(mean(Bias_squared)), Mean_Bias = mean(Bias), `90%_within` = quantile(Abs_Bias, probs = 0.9), Within_Interval = sum(Within_CI)/dplyr::n(), .groups = "drop_last") 
  
  return(list(estimates = estimates_out, summary_stats = summary_stats))
    
}