#' @title Summarize `rubias` Baseline Evaluation Tests
#'
#' @description
#' Thus function is a wrapper for [GCLr::custom_comb_rubias_output()] to summarize \pkg{rubias} output from baseline 
#' evaluation tests after running [GCLr::run_rubias_base_eval()].
#'
#' @param mixvec A character vector of test mixture names.
#' @param sample_sizes A tibble produced by [GCLr::base_eval_sample_sizes()] containing the following variables: 
#' `test_group`, `scenario`, `repunit`, and `samps`.
#' @param method A character vector of length 1 indicating the \pkg{rubias} output to summarize: 
#' "MCMC" (summarize MCMC output), "PB" (summarize bias corrected output), "both" (summarize both outputs); (default = "MCMC)
#' @param group_names An optional character vector of group names, used to sort `repunit` as a factor, passes through 
#' to [GCLr::custom_comb_rubias_output()] (default = `NULL`). If `NULL`, `group_names` comes from `repunit_trace.csv` output files.
#' @param group_names_new An optional character vector of new `group_names`, used to roll up groups from fine-scale `repunit` to broad-scale 
#' `repunit` for bias correction, passes through to [GCLr::custom_comb_rubias_output()] (default = `NULL`).
#' @param groupvec An optional numeric vector indicating the group affiliation of each pop in the baseline `sillyvec`, used if 
#' accessing `collection_trace.csv` output files to re-summarize to new groups (default = `NULL`). If `NULL`, `group_names` comes 
#' from `repunit_trace.csv` output files.
#' @param groupvec_new An optional numeric vector indicating the new broad-scale group affiliation of each fine-scale group, 
#' used if accessing `repunit_trace.csv` output files to re-summarize fine-scale groups to broad-scale groups with bias correction 
#' (`method = "PB"`; default = `NULL`). If `NULL`, `group_names` comes from `repunit_trace.csv` output files.
#' @param path A character vector of where to find output from each mixture as a .csv (created by [GCLr::run_rubias_base_eval()]; 
#' default is "rubias/output").
#' @param alpha A numeric vector of length 1 specifying credibility intervals (default is 0.1, which gives 90% CIs (i.e., 5% and 95%)).
#' @param burn_in A numeric vector of length 1 specifying how many sweeps were used for burn_in in [GCLr::run_rubias_base_eval()] 
#' (default = 5000).
#' @param threshold A numeric vector of length 1 specifying how low stock comp is before assume 0, used for `P=0` calculation, (default = 5e-7).
#' @param ncores An optional numeric value for the number of cores to use in a \pkg{foreach} `%dopar%` loop (default = 4). 
#' If the number of cores exceeds the number on your device ([parallel::detectCores()]), then all cores will be used. Note that `ncores` is 
#' only used if `method = "both"`.
#'
#' @returns A list with the following 2 components:
#'     \itemize{
#'       \item \code{estimates}: a tibble with 11 columns containing all group-level stock proportion output:
#'         \itemize{
#'           \item \code{test_group}: `repunit` being tested
#'           \item \code{scenario}: stock proportion of `repunit` being tested
#'           \item \code{repunit}: group name
#'           \item \code{mean}: mean posterior of stock proportion
#'           \item \code{sd}: sd of posterior of stock proportion
#'           \item \code{lo5CI}: lower 5% CI from posterior of stock proportion
#'           \item \code{hi95CI}: upper 95% CI from posterior of stock proportion 
#'           \item \code{P=0}: proportion of posterior with stock proportion < `threshold` (i.e., 0)
#'           \item \code{true_proportion}: the true stock proportion in  the scenario
#'           \item \code{total_samples}: mixture scenario sample size
#'           \item \code{method}: `rubias` method ("MCMC" or "PB")
#'         }
#'       \item \code{summary_stats}: a tibble with 6 columns containing mixture summary statistics for each group tested:
#'         \itemize{
#'           \item \code{method}: `rubias` method ("MCMC" or "PB")
#'           \item \code{test_group}: `repunit` being tested
#'           \item \code{RMSE}: root mean squared error of mean stock proportions
#'           \item \code{Mean_Bias}: mean bias (estimate - true) among all test mixtures for a given `test_group`
#'           \item \code{90_within}: the abs(90% CI) of the residuals (i.e., 90% of point estimates were within d distance of the true proportion)
#'           \item \code{Within_Interval}: the proportion of tests where the CI contained the true proportion
#'         }
#'     }
#'     
#' @details
#' Thus function is the final step in baseline evaluations tests using \pkg{rubias}. The normal workflow uses: [GCLr::base_eval_sample_sizes()] 
#' to determine leave-one-out sample sizes for tests, then [GCLr::create_rubias_base_eval()] generates the necessary \pkg{rubias} .csv files, which 
#' are then analyzed/run in parallel by [GCLr::run_rubias_base_eval()] to create output .csv files, and finally summarized by [GCLr::summarize_rubias_base_eval()].
#' 
#' @seealso 
#' [GCLr::custom_comb_rubias_output()]
#' [GCLr::base_eval_sample_sizes()]
#' [GCLr::create_rubias_base_eval()]
#' [GCLr::run_rubias_base_eval()]
#' [rubias::rubias()]
#' 
#' @examples
#' \dontrun{
#' base_eval_out <- GCLr::summarize_rubias_base_eval()
#' }
#' 
#' @export
summarize_rubias_base_eval <-
  function(mixvec,
           sample_sizes,
           method = c("MCMC", "PB", "both")[1],
           group_names = NULL,
           group_names_new = NULL,
           groupvec = NULL,
           groupvec_new = NULL,
           path = "rubias/output",
           alpha = 0.1,
           burn_in = 5000,
           threshold = 5e-7,
           ncores = 4) {

  start_time <- Sys.time()

  # Get summary estimates
  ## MCMC
  if (method == "MCMC") {
    estimates <-
      GCLr::custom_comb_rubias_output(
        rubias_output = NULL,
        mixvec = mixvec,
        group_names = group_names,
        group_names_new = group_names_new,
        groupvec = groupvec,
        groupvec_new = groupvec_new,
        path = path,
        alpha = alpha,
        burn_in = burn_in,
        bias_corr = FALSE,
        threshold = threshold,
        plot_trace = FALSE,
        ncores = ncores
      ) %>%
      dplyr::mutate(method = method)
    
  }
  
  ## PB (bias correction)
  if (method == "PB") {
    estimates <-
      GCLr::custom_comb_rubias_output(
        rubias_output = NULL,
        mixvec = mixvec,
        group_names = group_names,
        group_names_new = group_names_new,
        groupvec = groupvec,
        groupvec_new = groupvec_new,
        path = path,
        alpha = alpha,
        burn_in = burn_in,
        bias_corr = TRUE,
        threshold = threshold,
        plot_trace = FALSE,
        ncores = ncores
      ) %>%
      dplyr::mutate(method = method)
    
  }
  
  ## MCMC & PB
  if (method == "both") {
    
    if(ncores > parallel::detectCores()) {ncores = parallel::detectCores()}
    
    cl <- parallel::makePSOCKcluster(ncores)
    
    doParallel::registerDoParallel(cl, cores = ncores)
    
    parallel::clusterExport(cl = cl, varlist = "GCLr::custom_comb_rubias_output") # For some reason custom_comb_rubias_output() couldn't be found within the foreach loop, this solves that issue.
    
    `%dopar%` <- foreach::`%dopar%`
    
    estimates <-
      foreach::foreach(
        bias_corr = c(TRUE, FALSE),
        .packages = "tidyverse",
        .export = "GCLr:custom_comb_rubias_output"
      ) %dopar% {
        GCLr:custom_comb_rubias_output(
          rubias_output = NULL,
          mixvec = mixvec,
          group_names = group_names,
          group_names_new = group_names_new,
          groupvec = groupvec,
          groupvec_new = groupvec_new,
          path = path,
          alpha = alpha,
          burn_in = burn_in,
          bias_corr = bias_corr,
          threshold = threshold,
          plot_trace = FALSE,
          ncores = ncores
        ) %>%
          dplyr::mutate(method = if (bias_corr) {
            "PB"
          } else{
            "MCMC"
          })
        
      } %>% dplyr::bind_rows()
    
    parallel::stopCluster(cl)
    
  }
  
  # Create estimates output summary
  estimates_out <- estimates %>%
    tidyr::separate(
      mixture_collection,
      into = c("test_group", "scenario"),
      sep = "\\_(?=[^\\_]+$)",
      remove = TRUE
    ) %>% #The regular expression means the last instance of "_", just in case someone has a group name with an underscore.
    dplyr::mutate(scenario = as.numeric(scenario)) %>%
    dplyr::left_join(
      sample_sizes,
      by = c(
        "repunit" = "repunit",
        "scenario" = "scenario",
        "test_group" = "test_group"
      )
    ) %>%
    dplyr::group_by(test_group, scenario, method) %>%
    dplyr::mutate(total_samps = sum(samps)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      repunit = factor(repunit, levels = unique(repunit)),
      true_proportion = samps / total_samps,
      lo5CI = `5%`,
      hi95CI = `95%`
    ) %>%
    dplyr::select(
      test_group,
      scenario,
      repunit,
      mean,
      sd,
      lo5CI,
      hi95CI,
      `P=0`,
      true_proportion,
      total_samps,
      method
    )
  
  # Calculate the baseline eval summary statistics RMSE =  root mean squared error, 'Within_10%' =  proportion of estimates within 10% of true
  summary_stats <- estimates_out %>%
    dplyr::filter(test_group == repunit) %>%
    dplyr::group_by(method, test_group) %>%
    dplyr::mutate(
      Bias = mean - true_proportion,
      Abs_Bias = abs(mean - true_proportion),
      Bias_squared = (mean - true_proportion) ^ 2,
      CI_width = hi95CI - lo5CI,
      lo5CI_within = true_proportion - lo5CI,
      hi95CI_within = hi95CI - true_proportion,
      Within_CI = dplyr::case_when(
        true_proportion >= lo5CI &
          true_proportion <= hi95CI ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    dplyr::summarise(
      RMSE = sqrt(mean(Bias_squared)),
      Mean_Bias = mean(Bias),
      `90%_within` = quantile(Abs_Bias, probs = 0.9),
      Within_Interval = sum(Within_CI) / dplyr::n(),
      .groups = "drop_last"
    )
  
  return(list(estimates = estimates_out, summary_stats = summary_stats))
    
}