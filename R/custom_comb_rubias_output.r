#' @title Summarize `rubias` Output
#'
#' @description
#' A short description...
#'
#' @param rubias_output 
#' @param mixvec 
#' @param group_names 
#' @param group_names_new 
#' @param groupvec 
#' @param groupvec_new 
#' @param path 
#' @param alpha 
#' @param burn_in 
#' @param bias_corr 
#' @param threshold 
#' @param plot_trace 
#' @param ncores 
#'
#' @returns A tibble with 8 columns:
#'     \itemize{
#'       \item \code{mixture_collection}: factor of mixtures (only a factor for ordering, plotting purposes)
#'       \item \code{repunit}: factor of reporting groups (only a factor for ordering, plotting purposes)
#'       \item \code{mean}: mean posterior of stock proportion
#'       \item \code{sd}: sd of posterior of stock proportion
#'       \item \code{medoam}: median posterior of stock proportion
#'       \item \code{sd}: sd of posterior of stock proportion
#'       \item \code{loCI}: lower 5% CI from posterior of stock proportion
#'       \item \code{hiCI}: upper 95% CI from posterior of stock proportion 
#'       \item \code{P=0}: proportion of posterior with stock proportion < `threshold` (i.e., 0)
#'       }
#'
#' @examples
#' \dontrun{
#' rubias_out <- GCLr::custom_comb_rubias_output()
#' }
#' 
#' @export
custom_comb_rubias_output <-
  function(rubias_output = NULL,
           mixvec = NULL,
           group_names = NULL,
           group_names_new = NULL,
           groupvec = NULL,
           groupvec_new = NULL,
           path = "rubias/output",
           alpha = 0.1,
           burn_in = 5000,
           bias_corr = FALSE,
           threshold = 5e-7,
           plot_trace = TRUE,
           ncores = 4) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function computes summary statistics from `rubias` output, similar to `CustomCombineBAYESOutput.
  # However, output is a tibble with `mixture_collection` as a column, instead of each mixture as its own list.
  # It can take either the `rubias_output` list object from `run_rubias_mix` or `infer_mixture`,
  # OR it can read in the .csv files created by `run_rubias_mix`.
  #
  # NOTE: Currently this function only allows bias correction for the reporting groups run in the mixture
  # It can not do bias correction for different baseline groupvecs, because current `rubias` output only
  # gives the bias corrected means for each `mixture_collection` and `repunit` (i.e. `rho`, not `pi`)
  #
  # UPDATE: This function CAN do bias correction if you are rolling up groups from fine-scale to broad-scale.
  # To use this functionality, specify `group_names` as the original, fine-scale groups, `groupvec_new` as 
  # the groupvec to go from fine-scale to broad-scale groups, and `group_names_new` as the broad-scale groups
  # (i.e. length(groupvec_new) = length(group_names), and max(groupvec_new) == length(group_names_new))
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   rubias_output - output list object from `run_rubias_mix` or `infer_mixture`
  #   mixvec - character vector of mixture sillys, used to read in output .csv files if `rubias_output = NULL`
  #   group_names - character vector of group_names, used to sort repunit as a factor, can get from .csv
  #   group_names_new - character vector of new group_names, used to roll up groups from fine-scale to broad-scale for bias correction
  #                     if specified, `groupvec_new` must be = length(group_names_new)
  #   groupvec - numeric vector indicating the group affiliation of each pop in sillyvec, used if resumarizing to new groups
  #   groupvec_new - a numeric vector indicating the new group affiliation of each group, used if resumarizing fine-scale groups to broad-scale groups with bias correction
  #   path - character vector of where to find output from each mixture as a .csv (created by `run_rubias_mix`)
  #   alpha - numeric constant specifying credibility intervals, default is 0.1, which gives 90% CIs (i.e. 5% and 95%)
  #   burn_in - numeric constant specifying how many sweeps were used for burn_in in `run_rubias_mix` or `infer_mixture`
  #   bias_corr - logical switch indicating whether you want bias corrected values from `method = "PB"` or not, 
  #               currently can NOT do bias correction if not using the same repunits that were run in the mixture
  #   threshold - numeric constant specifying how low stock comp is before assume 0, used for `P=0` calculation, default is from BAYES
  #   plot_trace - logical switch, when on will create a trace plot for each mixture and repunit (reporting group)
  #   ncores - a numeric vector of length one indicating the number of cores to use (ncores is only used when is.null(rubias_output) == TRUE)
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a tibble with 8 fields for each mixture and repunit (reporting group)
  #   mixture_collection - factor of mixtures (only a factor for ordering, plotting purposes)
  #   repunit - factor of reporting groups (only a factor for ordering, plotting purposes)
  #   mean - mean stock composition
  #   sd - standard deviation
  #   median - median stock composition
  #   loCI - lower bound of credibility interval
  #   hiCI - upper bound of credibility interval
  #   P=0 - the proportion of the stock comp distribution that was below `threshold` (i.e posterior probability that stock comp = 0)
  #
  #   Also returns a trace plot for each mixture and repunit if `plot_trace = TRUE`
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load("V:/Analysis/1_SEAK/Sockeye/Mixture/Lynn Canal Inseason/2018/OLD rubias/output/test/custom_comb_rubias_output_test.RData")
  # lynncanal_2015_SW26_27.sum <- custom_comb_rubias_output(rubias_output = lynncanal_test.out, group_names = LynnCanal_groups7, bias_corr = TRUE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #~~~~~~~~~~~~~~~~
  ## Error catching
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  if(is.null(rubias_output) & !dir.exists(path)) {
    stop("`path` does not exist in your working directory, hoser!!!\nEither specify `rubias_output` or provide a valid `path` with rubias output .csv files for `mixvec`.")
  }
  if(is.null(rubias_output) & is.null(mixvec)) {
    stop("Need to provide either `rubias_output` tibble to summarize or `mixvec` and `path` so that rubias output can be read, hoser!!!")
  }
  if(!is.null(mixvec) & !all(sapply(mixvec, function(mixture) {any(grepl(pattern = mixture, x = list.files(path = path, pattern = ".csv")))} ))) {
    stop("Not all mixtures in `mixvec` have .csv output in `path`, hoser!!!")
  }
  if(!is.null(groupvec) & bias_corr) {
    stop("Can not perform bias correction if you are changing the groupvec (pop to group) from what was originally run.\n  Unfortunately since `rubias` only outputs bias corrected means for each `repunit`, we can't compute\n  bias corrected summary statistics on a new `groupvec`.")
  }
  if(!is.null(groupvec_new) & is.null(group_names_new)) {
    stop("Need to provide `group_names_new` if introducing `groupvec_new`, hoser!!!")
  }
  if(!is.null(groupvec_new)) {
    if(length(groupvec_new) != length(group_names) | max(groupvec_new) != length(group_names_new))
      stop("If specifying `groupvec_new`, you must be rolling up from fine-scale groups to broad-scale groups (i.e. length(groupvec_new) = length(group_names), and max(groupvec_new) == length(group_names_new))")
  }
  if(!is.null(groupvec) & is.null(group_names)) {
    stop("Need to provide `group_names` if introducing a new `groupvec`, hoser!!!")
  }
  
  `%dopar%` <- foreach::`%dopar%`
  
  #~~~~~~~~~~~~~~~~
  ## If no rubias_output, make from .csv files
  if(is.null(rubias_output)) {
    message("Summarizing results from rubias output .csv files.")
    
    # Build repunit_trace from "repunit_trace.csv" if `groupvec` is NULL
    if(is.null(groupvec)) {
      if(!all(file.exists(paste0(path, "/", mixvec, "_repunit_trace.csv")))) {
        stop("Not all mixtures in `mixvec` have a `repunit_trace.csv` file in `path`, hoser!!!")
      }  # make sure output files exist
      message("    Building trace output from `repunit_trace.csv` files.")
      
      cl <- parallel::makePSOCKcluster(ncores)
      
      doParallel::registerDoParallel(cl, cores = ncores)  
      
      repunit_trace <- foreach::foreach(mixture = mixvec, .packages = c("tidyverse")) %dopar% {
          
        repunit_trace_mix <- suppressMessages(readr::read_csv(file = paste0(path, "/", mixture, "_repunit_trace.csv")))
        repunit_trace_mix <- repunit_trace_mix %>% 
          tidyr::gather(repunit, rho, -sweep) %>%  # wide to tall
          dplyr::mutate(mixture_collection = mixture) %>%  # add mixture_collection
          dplyr::arrange(mixture_collection, sweep, repunit) %>%  # sort by mixture_collection, sweep, repunit
          dplyr::select(mixture_collection, sweep, repunit, rho)  # reorder columns
        
      } %>% dplyr::bind_rows()  # build output from "repunit_trace.csv"
      
      parallel::stopCluster(cl)
      
      if(is.null(group_names)) {
        group_names <- colnames(suppressMessages(readr::read_csv(file = paste0(path, "/", mixvec[1], "_repunit_trace.csv"))))[-1]
      }  # assign `group_names` from "repunit_trace.csv", if NULL
      
    } else {  # groupvec
      
      if(!all(file.exists(paste0(path, "/", mixvec, "_collection_trace.csv")))) {
        stop("Not all mixtures in `mixvec` have a `collection_trace.csv` file in `path`, hoser!!!")
      }  # make sure output files exist
      message("    Building trace output from `collection_trace.csv` files, `groupvec`, and `group_names`.")
      
      cl <- parallel::makePSOCKcluster(ncores)
      
      doParallel::registerDoParallel(cl, cores = ncores)  
      
      collection_trace <- foreach::foreach(mixture = mixvec, .packages = c("tidyverse")) %dopar% {
        
        collection_trace_mix <- suppressMessages(readr::read_csv(file = paste0(path, "/", mixture, "_collection_trace.csv")))
        collection_trace_mix <- collection_trace_mix %>% 
          tidyr::gather(collection, pi, -sweep) %>% 
          dplyr::mutate(mixture_collection = mixture)
        
      } %>% dplyr::bind_rows()  # build output from "collection_trace.csv"
      
      parallel::stopCluster(cl)
      
      base_collections <- unique(collection_trace$collection)  # baseline collections are the same order as in rubias output
      repunit_new.df <- tibble::tibble(collection = base_collections,
                                       repunit_new = group_names[groupvec])  # tibble of new repunit from groupvec
      repunit_trace <- collection_trace %>% 
        dplyr::left_join(repunit_new.df, by = "collection") %>%  # join with new repunit
        dplyr::rename(repunit = repunit_new) %>%  # rename new repunit
        dplyr::group_by(mixture_collection, sweep, repunit) %>%  # group and order
        dplyr::summarise(rho = sum(pi), .groups = "drop_last")   # summarise pi (collection) to rho (repunit)
        
    }  # build repunit_trace from "collection_trace.csv", `groupvec`, and `group_names`
    
    # Build `bootstrapped_proportions` if `bias_corr = TRUE`
    if(bias_corr) {
      if(!all(file.exists(paste0(path, "/", mixvec, "_bias_corr.csv")))) {
        stop("Not all mixtures in `mixvec` have a `bias_corr.csv` file in `path`, hoser!!!")
      }  # make sure output files exist
      message("    Building bias correction output from `bias_corr.csv` files.")

      cl <- parallel::makePSOCKcluster(ncores)
      
      doParallel::registerDoParallel(cl, cores = ncores)  
      
      bootstrapped_proportions <- foreach::foreach(mixture = mixvec, .packages = c("tidyverse")) %dopar% {
          
        bias_corr_mix <- suppressMessages(readr::read_csv(file = paste0(path, "/", mixture, "_bias_corr.csv")))
        bias_corr_mix <- bias_corr_mix %>% 
          dplyr::mutate(mixture_collection = mixture)
        
      } %>%  dplyr::bind_rows()  # build bootstrapped_proportions from "bias_corr.csv" files
      
      parallel::stopCluster(cl)
      
    }  # bias_corr
    
  }  # build rubias_output from .csv files, ignore "indiv_posteriors"
  
  #~~~~~~~~~~~~~~~~
  ## If rubias_output, create `repunit_trace` and `bootstrapped_proportions`
  if(!is.null(rubias_output)) {
    
    message("Summarizing results from `rubias output`.")
    
    if(is.null(mixvec)) {
      
      mixvec <- unique(rubias_output$mix_prop_traces$mixture_collection)
      
    }  # used to order mixtures as factor
    
    # Build repunit_trace from `rubias_output` if `groupvec` is NULL
    if(is.null(groupvec)) {
      
      repunit_trace <- rubias_output$mix_prop_traces %>% 
        dplyr::group_by(mixture_collection, sweep, repunit) %>%  # group to summarize across collections
        dplyr::summarise(rho = sum(pi), .groups = "drop_last" ) # summarize collections to repunits
        
      
      if(is.null(group_names)) {group_names <- unique(repunit_trace$repunit) }  # assign `group_names` from `rubias_output`, if NULL, order may be wrong
      
    } else {
      
      base_collections <- unique(rubias_output$mix_prop_traces$collection)  # baseline collections are the same order as in rubias output
      repunit_new.df <- tibble::tibble(collection = base_collections,
                                       repunit_new = group_names[groupvec])  # tibble of new repunit from groupvec
      repunit_trace <- rubias_output$mix_prop_traces %>% 
        dplyr::left_join(repunit_new.df, by = "collection") %>%  # join with new repunit
        dplyr::select(-repunit) %>%  # drop old repunit
        dplyr::rename(repunit = repunit_new) %>%  # rename new repunit
        dplyr::group_by(mixture_collection, sweep, repunit) %>%  # group and order
        dplyr::summarise(rho = sum(pi), .groups = "drop_last")   # summarise pi (collection) to rho (repunit)
        
    }  # build repunit_trace from `rubias_output`, `groupvec`, and `group_names`
    
    # Build `bootstrapped_proportions` if `bias_corr = TRUE`
    if(bias_corr) {
      message("    Building bias correction output from `rubias_output`.")
      bootstrapped_proportions <- rubias_output$bootstrapped_proportions
    }  # bias_corr
    
  }  # rubias_output
  
  #~~~~~~~~~~~~~~~~
  ## Calculate `d_rho` for bias correction if specified
  if(bias_corr) {
    if(nrow(bootstrapped_proportions) == 0) {stop("There is no bias corrected output, hoser!!!")}
    
    mixing_proportions_rho <- repunit_trace %>% 
      dplyr::filter(sweep >= burn_in) %>%   # remove burn_in
      dplyr::group_by(mixture_collection, repunit) %>%  # group by mixture and repunit across sweeps
      dplyr::summarise(rho = mean(rho), .groups = "drop_last") 
    
    d_rho <- mixing_proportions_rho %>% 
      dplyr::left_join(bootstrapped_proportions, by = c("mixture_collection", "repunit")) %>%  # join with `bootstrapped_proportions`
      dplyr::mutate(d_rho = rho - bs_corrected_repunit_ppn) %>%  # calculate d_rho
      dplyr::select(mixture_collection, repunit, d_rho)  # drop other variables
  }
  
  #~~~~~~~~~~~~~~~~
  ## Apply bias correction if `d_rho` exists
  if(exists("d_rho")) {
    repunit_trace <- repunit_trace %>% 
      dplyr::left_join(d_rho, by = c("mixture_collection", "repunit")) %>%  # join trace with d_rho
      dplyr::mutate(rho = rho - d_rho) %>%  # subtract d_rho
      dplyr::select(-d_rho)
  }
  
  #~~~~~~~~~~~~~~~~
  ## Roll up to broad-scale groups if `groupvec_new` specified
  if(!is.null(groupvec_new)) {
    level_key <- sapply(group_names, function(grp) {
      i = which(group_names == grp)
      group_names_new[groupvec_new[i]]
    }, simplify = FALSE )  # set up level_key to use with recode to roll up groups
    repunit_trace <- repunit_trace %>% 
      dplyr::mutate(repunit = recode(repunit, !!!level_key)) %>% 
      dplyr::group_by(mixture_collection, sweep, repunit) %>% 
      dplyr::summarise(rho = sum(rho), .groups = "drop_last") 
    grp_names <- group_names_new  # for factoring repunit
  } else {
    grp_names <- group_names  # for factoring repunit
  } 
  
  #~~~~~~~~~~~~~~~~
  ## Plot repunit trace
  if(plot_trace) {
    trace_plot <- repunit_trace %>% 
      dplyr::mutate(mixture_collection = factor(x = mixture_collection, levels = mixvec)) %>%  # order mixture_collection
      dplyr::mutate(repunit = factor(x = repunit, levels = grp_names)) %>%  # order repunit
      ggplot2::ggplot(aes(x = sweep, y = rho, colour = repunit)) +
      ggplot2::geom_line() +
      ggplot2::ylim(0, 1) +
      ggplot2::geom_vline(xintercept = burn_in) +
      # ggplot2::annotate(geom = "text", x = burn_in / 2, y = 0.9, label = "Burn-in") +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::facet_grid(repunit ~ mixture_collection)
  }

  #~~~~~~~~~~~~~~~~
  ## Summary statistics
  loCI = alpha / 2
  hiCI = 1 - (alpha / 2)
  
  out_sum <- repunit_trace %>% 
    dplyr::filter(sweep >= burn_in) %>%  # remove burn_in
    dplyr::mutate(mixture_collection = factor(x = mixture_collection, levels = mixvec)) %>%  # order mixture_collection
    dplyr::mutate(repunit = factor(x = repunit, levels = grp_names)) %>%  # order repunit
    dplyr::group_by(mixture_collection, repunit) %>%  # group by mixture and repunit across sweeps
    dplyr::summarise(mean = mean(rho),
                     sd = sd(rho),
                     median = median(rho),
                     loCI = quantile(rho, probs = loCI),
                     hiCI = quantile(rho, probs = hiCI),
                     `P=0` = sum(rho < threshold) / length(rho), .groups = "drop_last") %>%    # summary statistics to return
    dplyr::mutate(loCI = replace(loCI, which(loCI < 0), 0),
                  hiCI = replace(hiCI, which(hiCI < 0), 0),
                  median = replace(median, which(median < 0), 0),
                  mean = replace(mean, which(mean < 0), 0)) %>% 
    dplyr::mutate(loCI = replace(loCI, which(loCI > 1), 1),
                  hiCI = replace(hiCI, which(hiCI > 1), 1),
                  median = replace(median, which(median > 1), 1),
                  mean = replace(mean, which(mean > 1), 1)) %>% 
    magrittr::set_colnames(c("mixture_collection", "repunit", "mean", "sd", "median", 
                             paste0(loCI * 100, "%"), paste0(hiCI * 100, "%"), "P=0"))
  
  if(exists("trace_plot")) {suppressWarnings(print(trace_plot))}  # plot the trace for each mixture, repunit
  return(out_sum)
}  # end function