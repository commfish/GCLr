#' @title Summarize `rubias` Baseline Evaluation Tests
#'
#' @description
#' Thus function summarizes \pkg{rubias} output from baseline 
#' evaluation tests after running [GCLr::run_rubias_base_eval()].
#'
#' @param mixvec A character vector of test mixture names.
#' @param sample_sizes A tibble produced by [GCLr::base_eval_sample_sizes()] containing the following variables: 
#' `test_group`, `scenario`, `repunit`, and `samps`.
#' @param method A character vector of length 1 indicating the \pkg{rubias} output to summarize: 
#' "MCMC" (summarize MCMC output), "PB" (summarize bias corrected output), "both" (summarize both outputs); (default = "MCMC")
#' @param path A character vector of where to find output from each mixture as a .csv or .fst file (created by [GCLr::run_rubias_base_eval()]; 
#' default is "rubias/output").
#' @param alpha A numeric vector of length 1 specifying credibility intervals (default is 0.1, which gives 90% CIs (i.e., 5% and 95%)).
#' @param burn_in A numeric vector of length 1 specifying how many sweeps were used for burn_in in [GCLr::run_rubias_base_eval()] 
#' (default = 5000).
#' @param threshold A numeric vector of length 1 specifying how low stock comp is before assume 0, used for `P=0` calculation, (default = 5e-7).
#' @param ncores An optional numeric value for the number of cores to use in a \pkg{foreach} `%dopar%` loop (default = 4). 
#' If the number of cores exceeds the number on your device ([parallel::detectCores()]), then all cores will be used. Note that `ncores` is 
#' only used if `method = "both"`.
#' @param file_type the file type of the mixture output files .fst (default and faster) or .csv files.
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
#' to determine test mixture sample sizes, then [GCLr::create_rubias_base_eval()] generates the necessary \pkg{rubias} files, which 
#' are then analyzed/run in parallel by [GCLr::run_rubias_base_eval()] to create output files, and finally summarized by [GCLr::summarize_rubias_base_eval()].
#' This function has the option of summarizing fst or csv mixture output files. fst files are compressed so they read into R faster, which speeds up the mixture summary process.
#' 
#' @seealso 
#' [GCLr::base_eval_sample_sizes()]
#' [GCLr::create_rubias_base_eval()]
#' [GCLr::run_rubias_base_eval()]
#' [rubias::rubias()]
#' 
#' @examples
#' \dontrun{
#'  
#' path <- "V:/Analysis/5_Coastwide/Sockeye/IYS_2022_MSA/rubias/Evaluating IYS sockeye baseline groups/rubias/output"
#' 
#' sample_sizes <- readRDS("V:/Analysis/5_Coastwide/Sockeye/IYS_2022_MSA/rubias/Evaluating IYS sockeye baseline groups/output/sample_sizes.rds")
#' 
#' mixvec <- suppressWarnings(tibble::tibble(files = list.files(path, pattern = "_repunit_trace.csv")) %>% 
#'                              tidyr::separate(files, into = c("mixture", NA), sep = "_repunit_trace.csv") %>% 
#'                              dplyr::pull(mixture))
#' 
#' eval_output <- GCLr::summarize_rubias_base_eval(mixvec = mixvec, path = path, sample_sizes = sample_sizes, ncores = parallel::detectCores(), file_type = "csv")
#' 
#' }
#' 
#' @export
summarize_rubias_base_eval <- function(
                                       mixvec,
                                       sample_sizes,
                                       method = c("MCMC", "PB", "both")[1],
                                       path = "rubias/output",
                                       alpha = 0.1,
                                       burn_in = 5000,
                                       threshold = 5e-7,
                                       ncores = 4,
                                       file_type = c("fst", "csv")[1]) {
  start_time <- Sys.time()
  
  # Error catching ----
  if (ncores > parallel::detectCores()) {
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine"    )
  }
  
  if (length(list.files(path = path, pattern = paste0(".", file_type))) == 0) {
    stop(paste0("There are no files with `file_type` .", file_type, " output in `path`, `file_type` must be 'csv' or 'fst'"))
  }  # catch file_type errors before checking all `mixvec` since that step can be slow
  
  if (!is.null(mixvec) & !all(sapply(mixvec, function(mixture) {any(grepl(pattern = mixture, x = list.files(path = path, pattern = paste0(".", file_type))))} ))) {
    stop(paste0("Not all mixtures in `mixvec` have .", file_type, " output in `path`, hoser!!!"))
  }
  
      p_objects <-
        c(
          "mixvec",
          "sample_sizes",
          "method",
          "path",
          "alpha",
          "burn_in",
          "threshold",
          "ncores",
          "file_type"
        ) #Objects needed inside the parallel loop.
      
      # Get summary estimates
      ## MCMC
      if (method == "MCMC") {
        message(paste0("    Building MCMC output from ", file_type, " files."))
        
        cl <- parallel::makePSOCKcluster(ncores)
        
        parallel::clusterExport(cl = cl,
                                varlist = p_objects,
                                envir = environment())
        parallel::clusterEvalQ(cl = cl, library(tidyverse))
        
        estimates_out <-
          pbapply::pblapply(
            cl = cl,
            X = 1:length(mixvec),
            FUN = function(i) {
              mixture <- mixvec[i]
              
              if (file_type == "csv") {
                repunit_trace_mix <-
                  suppressMessages(readr::read_csv(file = paste0(
                    path, "/", mixture, "_repunit_trace.csv"
                  )))
                
              } else {
                repunit_trace_mix <-
                  suppressMessages(fst::read_fst(path = paste0(
                    path, "/", mixture, "_repunit_trace.fst"
                  )))
                
              }
              
              repunit_trace <- repunit_trace_mix %>%
                dplyr::select(-tidyselect::any_of("chain")) %>%  # old functions didn't include multichain
                tidyr::pivot_longer(-sweep,
                                    names_to = "repunit",
                                    values_to = "rho") %>%  # wide to tall
                dplyr::mutate(mixture_collection = mixture) %>%
                dplyr::arrange(mixture_collection, sweep, repunit) %>%
                dplyr::select(mixture_collection, sweep, repunit, rho)  # reorder columns
              
              grp_names <- repunit_trace$repunit %>% unique()
              
              #~~~~~~~~~~~~~~~~
              ## Summary statistics ----
              loCI = alpha / 2
              hiCI = 1 - (alpha / 2)
              
              estimates <- repunit_trace %>%
                dplyr::filter(sweep >= burn_in) %>%  # remove burn_in
                dplyr::mutate(
                  mixture_collection = factor(x = mixture_collection, levels = mixvec),
                  repunit = factor(x = repunit, levels = grp_names)
                ) %>%  # order mixture_collection, repunit
                dplyr::group_by(mixture_collection, repunit) %>%
                dplyr::summarise(
                  mean = mean(rho),
                  sd = sd(rho),
                  median = median(rho),
                  loCI = quantile(rho, probs = loCI),
                  hiCI = quantile(rho, probs = hiCI),
                  `P=0` = sum(rho < threshold) / length(rho),
                  .groups = "drop"
                ) %>%
                dplyr::mutate(
                  loCI = replace(loCI, which(loCI < 0), 0),
                  hiCI = replace(hiCI, which(hiCI < 0), 0),
                  median = replace(median, which(median < 0), 0),
                  mean = replace(mean, which(mean < 0), 0)
                ) %>%
                dplyr::mutate(
                  loCI = replace(loCI, which(loCI > 1), 1),
                  hiCI = replace(hiCI, which(hiCI > 1), 1),
                  median = replace(median, which(median > 1), 1),
                  mean = replace(mean, which(mean > 1), 1)
                ) %>%
                magrittr::set_colnames(
                  c(
                    "mixture_collection",
                    "repunit",
                    "mean",
                    "sd",
                    "median",
                    paste0(loCI * 100, "%"),
                    paste0(hiCI * 100, "%"),
                    "P=0",
                    "GR",
                    "n_eff"
                  )
                ) %>%
                dplyr::mutate(method = "MCMC") %>%
                dplyr::select(method, dplyr::everything())
              
              # Create estimates output summary
              estimates %>%
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
              
            }
          ) %>% dplyr::bind_rows()  # Get estimates from trace files
        
        parallel::stopCluster(cl)
        
      } #MCMC only
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      ## PB (bias correction)
      if (method == "PB") {
        if (!all(file.exists(paste0(
          path, "/", mixvec, "_bias_corr.", file_type
        )))) {
          stop(
            paste0(
              "Not all mixtures in `mixvec` have a `bias_corr.",
              file_type,
              "` file in `path`, hoser!!!"
            )
          )
        }  # make sure output files exist
        
        message(paste0("    Building PB output from ", file_type, " files."))
        
        cl <- parallel::makePSOCKcluster(ncores)
        
        parallel::clusterExport(cl = cl,
                                varlist = p_objects,
                                envir = environment())
        parallel::clusterEvalQ(cl = cl, library(tidyverse))
        
        estimates_out <-
          pbapply::pblapply(
            cl = cl,
            X = 1:length(mixvec),
            FUN = function(i) {
              mixture <- mixvec[i]
              
              if (file_type == "csv") {
                repunit_trace_mix <-
                  suppressMessages(readr::read_csv(file = paste0(
                    path, "/", mixture, "_repunit_trace.csv"
                  )))
                
              } else {
                repunit_trace_mix <-
                  suppressMessages(fst::read_fst(path = paste0(
                    path, "/", mixture, "_repunit_trace.fst"
                  )))
                
              }
              
              repunit_trace <- repunit_trace_mix %>%
                dplyr::select(-tidyselect::any_of("chain")) %>%  # old functions didn't include multichain
                tidyr::pivot_longer(-sweep,
                                    names_to = "repunit",
                                    values_to = "rho") %>%  # wide to tall
                dplyr::mutate(mixture_collection = mixture) %>%
                dplyr::arrange(mixture_collection, sweep, repunit) %>%
                dplyr::select(mixture_collection, sweep, repunit, rho)  # reorder columns
              
              # Get bias corrected estimates
              if (file_type == "csv") {
                bias_corr <-
                  suppressMessages(readr::read_csv(file = paste0(
                    path, "/", mixture, "_bias_corr.csv"
                  )))
                
              } else {
                bias_corr <-
                  suppressMessages(fst::read_fst(path = paste0(
                    path, "/", mixture, "_bias_corr.fst"
                  )))
                
              }
              
              bootstrapped_proportions <- bias_corr %>%
                dplyr::mutate(mixture_collection = mixture)
              
              ## Calculate `d_rho` for bias correction
              
              mixing_proportions_rho <- repunit_trace %>%
                dplyr::filter(sweep >= burn_in) %>%   # remove burn_in
                dplyr::group_by(mixture_collection, repunit) %>%  # group by mixture and repunit across sweeps
                dplyr::summarise(rho = mean(rho), .groups = "drop")
              
              d_rho <- mixing_proportions_rho %>%
                dplyr::left_join(bootstrapped_proportions,
                                 by = c("mixture_collection", "repunit")) %>%  # join with `bootstrapped_proportions`
                dplyr::mutate(d_rho = rho - bs_corrected_repunit_ppn) %>%  # calculate d_rho
                dplyr::select(mixture_collection, repunit, d_rho)  # drop other variables
              
              ## Bias correct the MCMC output
              repunit_trace <- repunit_trace %>%
                dplyr::left_join(d_rho, by = c("mixture_collection", "repunit")) %>%  # join trace with d_rho
                dplyr::mutate(rho = rho - d_rho) %>%  # subtract d_rho
                dplyr::select(-d_rho)
              
              grp_names <- repunit_trace$repunit %>% unique()
              
              #~~~~~~~~~~~~~~~~
              ## Summary statistics ----
              loCI = alpha / 2
              hiCI = 1 - (alpha / 2)
            
              estimates <- repunit_trace %>%
                dplyr::filter(sweep >= burn_in) %>%  # remove burn_in
                dplyr::mutate(
                  mixture_collection = factor(x = mixture_collection, levels = mixvec),
                  repunit = factor(x = repunit, levels = grp_names)
                ) %>%  # order mixture_collection, repunit
                dplyr::group_by(mixture_collection, repunit) %>%
                dplyr::summarise(
                  mean = mean(rho),
                  sd = sd(rho),
                  median = median(rho),
                  loCI = quantile(rho, probs = loCI),
                  hiCI = quantile(rho, probs = hiCI),
                  `P=0` = sum(rho < threshold) / length(rho),
                  .groups = "drop"
                ) %>%
                dplyr::mutate(
                  loCI = replace(loCI, which(loCI < 0), 0),
                  hiCI = replace(hiCI, which(hiCI < 0), 0),
                  median = replace(median, which(median < 0), 0),
                  mean = replace(mean, which(mean < 0), 0)
                ) %>%
                dplyr::mutate(
                  loCI = replace(loCI, which(loCI > 1), 1),
                  hiCI = replace(hiCI, which(hiCI > 1), 1),
                  median = replace(median, which(median > 1), 1),
                  mean = replace(mean, which(mean > 1), 1)
                ) %>%
                magrittr::set_colnames(
                  c(
                    "mixture_collection",
                    "repunit",
                    "mean",
                    "sd",
                    "median",
                    paste0(loCI * 100, "%"),
                    paste0(hiCI * 100, "%"),
                    "P=0"
                  )
                ) %>%
                dplyr::mutate(method = "PB") %>%
                dplyr::select(method, dplyr::everything())
              
              # Create estimates output summary
              estimates %>%
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
              
            }
          ) %>% dplyr::bind_rows()  # Get estimates from trace files
        
        parallel::stopCluster(cl)
        
      } #PB only
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      ## MCMC & PB
      if (method == "both") {
        if (!all(file.exists(paste0(
          path, "/", mixvec, "_bias_corr.", file_type
        )))) {
          stop(
            paste0(
              "Not all mixtures in `mixvec` have a `bias_corr.",
              file_type,
              "` file in `path`, hoser!!!"
            )
          )
        }  # make sure output files exist
        
        message(paste0("    Building MCMC and PB output from ", file_type, " files."))
        
        cl <- parallel::makePSOCKcluster(ncores)
        
        parallel::clusterExport(cl = cl,
                                varlist = p_objects,
                                envir = environment())
        parallel::clusterEvalQ(cl = cl, library(tidyverse))
        
        estimates_out <-
          pbapply::pblapply(
            cl = cl,
            X = 1:length(mixvec),
            FUN = function(i) {
              
              mixture <- mixvec[i]
              
              if (file_type == "csv") {
                repunit_trace_mix <-
                  suppressMessages(readr::read_csv(file = paste0(
                    path, "/", mixture, "_repunit_trace.csv"
                  )))
                
              } else{
                repunit_trace_mix <-
                  suppressMessages(fst::read_fst(path = paste0(
                    path, "/", mixture, "_repunit_trace.fst"
                  )))
                
              }
              
              repunit_trace <- repunit_trace_mix %>%
                dplyr::select(-tidyselect::any_of("chain")) %>%  # old functions didn't include multichain
                tidyr::pivot_longer(-sweep,
                                    names_to = "repunit",
                                    values_to = "rho") %>%  # wide to tall
                dplyr::mutate(mixture_collection = mixture) %>%
                dplyr::arrange(mixture_collection, sweep, repunit) %>%
                dplyr::select(mixture_collection, sweep, repunit, rho)  # reorder columns
              
              # Get bias corrected estimates
              if (file_type == "csv") {
                bias_corr <-
                  suppressMessages(readr::read_csv(file = paste0(
                    path, "/", mixture, "_bias_corr.csv"
                  )))
                
              } else {
                bias_corr <-
                  suppressMessages(fst::read_fst(path = paste0(
                    path, "/", mixture, "_bias_corr.fst"
                  )))
                
              }
              
              bootstrapped_proportions <- bias_corr %>%
                dplyr::mutate(mixture_collection = mixture)
              
              ## Calculate `d_rho` for bias correction
              
              mixing_proportions_rho <- repunit_trace %>%
                dplyr::filter(sweep >= burn_in) %>%   # remove burn_in
                dplyr::group_by(mixture_collection, repunit) %>%  # group by mixture and repunit across sweeps
                dplyr::summarise(rho = mean(rho), .groups = "drop")
              
              d_rho <- mixing_proportions_rho %>%
                dplyr::left_join(bootstrapped_proportions,
                                 by = c("mixture_collection", "repunit")) %>%  # join with `bootstrapped_proportions`
                dplyr::mutate(d_rho = rho - bs_corrected_repunit_ppn) %>%  # calculate d_rho
                dplyr::select(mixture_collection, repunit, d_rho)  # drop other variables
              
              ## Bias correct the MCMC output
              repunit_trace <- repunit_trace %>%
                dplyr::left_join(d_rho, by = c("mixture_collection", "repunit")) %>%  # join trace with d_rho
                dplyr::mutate(pb_rho = rho - d_rho) %>%  # subtract d_rho
                dplyr::select(-d_rho)
              
              grp_names <- repunit_trace$repunit %>% unique()
              
              #~~~~~~~~~~~~~~~~
              ## Summary statistics ----
              loCI = alpha / 2
              hiCI = 1 - (alpha / 2)
              
              estimates_mcmc <- repunit_trace %>%
                select(-pb_rho) %>%
                dplyr::filter(sweep >= burn_in) %>%  # remove burn_in
                dplyr::mutate(
                  mixture_collection = factor(x = mixture_collection, levels = mixvec),
                  repunit = factor(x = repunit, levels = grp_names)
                ) %>%  # order mixture_collection, repunit
                dplyr::group_by(mixture_collection, repunit) %>%
                dplyr::summarise(
                  mean = mean(rho),
                  sd = sd(rho),
                  median = median(rho),
                  loCI = quantile(rho, probs = loCI),
                  hiCI = quantile(rho, probs = hiCI),
                  `P=0` = sum(rho < threshold) / length(rho),
                  .groups = "drop"
                ) %>%
                dplyr::mutate(
                  loCI = replace(loCI, which(loCI < 0), 0),
                  hiCI = replace(hiCI, which(hiCI < 0), 0),
                  median = replace(median, which(median < 0), 0),
                  mean = replace(mean, which(mean < 0), 0)
                ) %>%
                dplyr::mutate(
                  loCI = replace(loCI, which(loCI > 1), 1),
                  hiCI = replace(hiCI, which(hiCI > 1), 1),
                  median = replace(median, which(median > 1), 1),
                  mean = replace(mean, which(mean > 1), 1)
                ) %>%
                magrittr::set_colnames(
                  c(
                    "mixture_collection",
                    "repunit",
                    "mean",
                    "sd",
                    "median",
                    paste0(loCI * 100, "%"),
                    paste0(hiCI * 100, "%"),
                    "P=0"
                  )
                ) %>%
                dplyr::mutate(method = "MCMC") %>%
                dplyr::select(method, dplyr::everything())
              
              
              # Bias corrected (PB)
              estimates_pb <- repunit_trace  %>%
                mutate(rho = pb_rho) %>%
                select(-pb_rho) %>%
                dplyr::filter(sweep >= burn_in) %>%  # remove burn_in
                dplyr::mutate(
                  mixture_collection = factor(x = mixture_collection, levels = mixvec),
                  repunit = factor(x = repunit, levels = grp_names)
                ) %>%  # order mixture_collection, repunit
                dplyr::group_by(mixture_collection, repunit) %>%
                dplyr::summarise(
                  mean = mean(rho),
                  sd = sd(rho),
                  median = median(rho),
                  loCI = quantile(rho, probs = loCI),
                  hiCI = quantile(rho, probs = hiCI),
                  `P=0` = sum(rho < threshold) / length(rho),
                  .groups = "drop"
                ) %>%
                dplyr::mutate(
                  loCI = replace(loCI, which(loCI < 0), 0),
                  hiCI = replace(hiCI, which(hiCI < 0), 0),
                  median = replace(median, which(median < 0), 0),
                  mean = replace(mean, which(mean < 0), 0)
                ) %>%
                dplyr::mutate(
                  loCI = replace(loCI, which(loCI > 1), 1),
                  hiCI = replace(hiCI, which(hiCI > 1), 1),
                  median = replace(median, which(median > 1), 1),
                  mean = replace(mean, which(mean > 1), 1)
                ) %>%
                magrittr::set_colnames(
                  c(
                    "mixture_collection",
                    "repunit",
                    "mean",
                    "sd",
                    "median",
                    paste0(loCI * 100, "%"),
                    paste0(hiCI * 100, "%"),
                    "P=0"
                  )
                ) %>%
                dplyr::mutate(method = "PB") %>%
                dplyr::select(method, dplyr::everything()) # PB estimates
              
              #Combine MCMC and bias corrected estimates
              estimates <- bind_rows(estimates_mcmc, estimates_pb)
              
              # Create estimates output summary
              estimates %>%
                dplyr::bind_rows() %>%
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
              
            }) %>% dplyr::bind_rows()
        
        parallel::stopCluster(cl)
        
      }# Both
      
      # Calculate the baseline eval summary statistics RMSE =  root mean squared error, 'Within_10%' =  proportion of estimates within 10% of true
      
      message("Calulating summary statistics")
      
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
      
      Sys.time() - start_time
      
      return(list(estimates = estimates_out, summary_stats = summary_stats))
      
}
