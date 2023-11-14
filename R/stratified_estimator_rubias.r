#' Summary statistics from a stratified estimate of rubias output
#'
#' This function computes summary statistics from a stratified estimate of \pkg{rubias} output,
#' as a column, instead of a single matrix. It can take either the `rubias_output` list object
#' from `run_rubias_mix()` or [rubias::infer_mixture()], OR it can read in the .csv files created by `run_rubias_mix()`.
#'
#' NOTE: Currently this function only allows bias correction for the reporting groups run in the mixture.
#' It cannot do bias correction for different baseline groupvecs, because current \pkg{rubias} output only
#' gives the bias-corrected means for each `mixture_collection` and `repunit` (i.e., `rho`, not `pi`).
#'
#' UPDATE:
#'  - This function CAN do bias correction if you are rolling up groups from fine-scale to broad-scale. To use this functionality, specify `group_names` as the original, fine-scale groups, `groupvec_new` as the groupvec to go from fine-scale to broad-scale groups, and `group_names_new` as the broad-scale groups (i.e., length(groupvec_new) = length(group_names), and max(groupvec_new) == length(group_names_new)).
#'  - 3/18/2022 corrected some errors in calculations (only when catch numbers had cv's). New output format now includes both catch numbers and stock props columns with '_backcalc' are the stock props back-calculated using simulated harvest stock props are cut off at 1.
#'  - Multi-chain compatible
#' 
#' @param rubias_output Output list object from `run_rubias_mix()` or [rubias::infer_mixture()].
#' @param mixvec Character vector of mixture sillies, used to read in output .csv files if `rubias_output = NULL`.
#' @param group_names Character vector of group_names, used to sort repunit as a factor, can get from .csv.
#' @param catchvec Numeric vector of harvest for each strata, must be in the same order as `mixvec`.
#' @param newname Character vector of length 1 specifying the name of the stratified estimate.
#' @param group_names_new Character vector of new group_names, used to roll up groups from fine-scale to broad-scale for bias correction.
#' @param groupvec Numeric vector indicating the group affiliation of each pop in sillyvec, used if resuming to new groups.
#' @param groupvec_new A numeric vector indicating the new group affiliation of each group, used if resuming fine-scale groups to broad-scale groups with bias correction.
#' @param path Character vector of where to find output from each mixture as a .csv (created by `run_rubias_mix()`).
#' @param alpha Numeric constant specifying credibility intervals, default is 0.1, which gives 90% CIs (i.e., 5% and 95%).
#' @param burn_in Numeric constant specifying how many sweeps were used for burn_in in `run_rubias_mix()` or [rubias::infer_mixture()].
#' @param bias_corr Logical switch indicating whether you want bias-corrected values from `method = "PB"` or not, 
#'                  currently can NOT do bias correction if not using the same repunits that were run in the mixture.
#' @param threshold Numeric constant specifying how low stock comp is before assuming 0, used for `P=0` calculation, default is from BAYES.
#' @param cv Numeric vector of harvest estimate coefficients of variation for each stratum, must be the same order as `mixvec`. 
#'           (Only for catch numbers with CVs).
#'
#' @return Returns a tibble with 13 fields (18 if with CVs) for each repunit (reporting group).
#'   - stratified_mixture: Character of stratified mixture.
#'   - repunit: Factor of reporting groups (only a factor for ordering, plotting purposes).
#'   - mean_harv: Mean stock catch number.
#'   - sd_harv: Standard deviation for catch.
#'   - median_harv: Median stock catch number.
#'   - `5%_harv` (depends on alpha level): Lower bound of credibility interval for catch number.
#'   - `95%_harv`: Upper bound of credibility interval for catch number.
#'   - mean: Mean stock composition.
#'   - sd: Standard deviation.
#'   - median: Median stock composition.
#'   - `5%`: Lower bound of credibility interval.
#'   - `95%`: Upper bound of credibility interval.
#'   - P=0: The proportion of the stock comp distribution that was below threshold (i.e., posterior probability that stock comp = 0).
#'   (For catch numbers with CVs):
#'   - mean_backcalc: mean_harv/sum(mean_harv); mean stock composition.
#'   - sd_backcalc: sd_harv/sum(mean_harv); standard deviation.
#'   - median_backcalc: median_harv/sum(mean_harv); median stock composition.
#'   - `5%_backcalc`: `5%_harv`/sum(mean_harv); lower bound of credibility interval.
#'   - `95%_backcalc`: `95%_harv`/sum(mean_harv); upper bound of credibility interval.
#'
#' @examples
#' path <- "V:/Analysis/1_SEAK/Sockeye/Mixture/D101to103 Seine/output_PB"
#' 
#' GCLr::stratified_estimator_rubias(path = path,
#'                             mixvec = c("D101__28293031", "D101__323334", "D101__3536"),
#'                             group_names = c("Alaska", "Nass", "Skeena", "Other", "McDonald", "Hugh Smith", "Klawock"),
#'                             catchvec = c(34860, 38169, 21562),
#'                             newname = "D101_21_stratified_NBDomestic",
#'                             group_names_new = c("Alaska", "Nass", "Skeena", "Other", "McDonald", "Hugh Smith"),
#'                             groupvec_new = c(1, 2, 3, 4, 5, 6, 1),
#'                             alpha = 0.1,
#'                             burn_in = 5000,
#'                             bias_corr = TRUE,
#'                             cv = c(0.2, 0.5, 0.4))
#' 
#' @export
stratified_estimator_rubias <-
  function(rubias_output = NULL,
           mixvec = NULL,
           group_names = NULL,
           catchvec, newname = NULL,
           group_names_new = NULL,
           groupvec = NULL,
           groupvec_new = NULL,
           path = "rubias/output",
           alpha = 0.1,
           burn_in = 5000,
           bias_corr = FALSE,
           threshold = 5e-7,
           cv = NULL,
           ncores = 4) {
  
  # Error catching ----
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  if(is.null(rubias_output) & !dir.exists(path)){
    
    stop("`path` does not exist in your working directory!!!\nEither specify `rubias_output` or provide a valid `path` with rubias output .csv files for `mixvec`.")
    
  }
  
  if(is.null(rubias_output) & is.null(mixvec)) {
    
    stop("Need to provide either `rubias_output` tibble to summarize or `mixvec` and `path` so that rubias output can be read!!!")
    
  }
  
  if(!is.null(mixvec) & is.null(rubias_output) & !all(sapply(mixvec, function(mixture) {any(grepl(pattern = mixture, x = list.files(path = path, pattern = ".csv")))} ))) {
    
    stop("Not all mixtures in `mixvec` have .csv output in `path`!!!")
    
  }
  
  if(!is.null(groupvec) & bias_corr) {
    
    stop("Can not perform bias correction if you are changing the groupvec (pop to group) from what was originally run.\n  Unfortunately since `rubias` only outputs bias corrected means for each `repunit`, we can't compute\n  bias corrected summary statistics on a new `groupvec`.")
    
  }
  
  if(!is.null(groupvec_new) & is.null(group_names_new)) {
    
    stop("Need to provide `group_names_new` if introducing `groupvec_new`!!!")
    
  }
  
  if(!is.null(groupvec_new)) {
    
    if(length(groupvec_new) != length(group_names) | max(groupvec_new) != length(group_names_new))
      
      stop("If specifying `groupvec_new`, you must be rolling up from fine-scale groups to broad-scale groups (i.e. length(groupvec_new) = length(group_names), and max(groupvec_new) == length(group_names_new))")
    
  }
  
  if(!is.null(groupvec) & is.null(group_names)) {
    
    stop("Need to provide `group_names` if introducing a new `groupvec`")
    
  }
  
  if(is.null(newname)) {
    
    stop("Need to provide a `newname` for this stratified estimate!!!")
    
  }
  
  `%dopar%` <- foreach::`%dopar%`
  
  # Summaries ----
  ## If no rubias_output, ----
  # make from .csv files
  if(is.null(rubias_output)) {
    
    message("Summarizing results from rubias output .csv files.")
    
    # Build repunit_trace from "repunit_trace.csv" if `groupvec` is NULL
    if(is.null(groupvec)) {
      
      if(!all(file.exists(paste0(path, "/", mixvec, "_repunit_trace.csv")))) {
        
        stop("Not all mixtures in `mixvec` have a `repunit_trace.csv` file in `path`!!!")
        
      }  # make sure output files exist
      
      message("    Building trace output from `repunit_trace.csv` files.")
      
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl, cores = ncores)  
      
      repunit_trace <-
        foreach::foreach(mixture = mixvec, .packages = c("tidyverse")) %dopar% {
          
          repunit_trace_mix <-
            suppressMessages(readr::read_csv(file = paste0(path, "/", mixture, "_repunit_trace.csv")))
          
          if (!"chain" %in% names(repunit_trace_mix)) {
            dplyr::mutate(repunit_trace_mix, chain = 1L)
          } # in case older files without multichain
          
          repunit_trace_mix <- repunit_trace_mix %>% 
            tidyr::pivot_longer(-c(sweep, chain), names_to = "repunit", values_to = "rho") %>%  # wide to tall
            dplyr::mutate(mixture_collection = mixture) %>%
            dplyr::arrange(mixture_collection, chain, sweep, repunit) %>%
            dplyr::select(mixture_collection, chain, sweep, repunit, rho)  # reorder columns
          
        } %>% dplyr::bind_rows()  # build output from "repunit_trace.csv"
      
      parallel::stopCluster(cl)
      
      if(is.null(group_names)) {
        
        group_names <-
          suppressMessages(readr::read_csv(file = paste0(path, "/", mixvec[1], "_repunit_trace.csv"))) %>%
          colnames() %>%
          {.[which(!. %in% c("sweep", "chain"))]}
        
      }  # assign `group_names` from "repunit_trace.csv", if NULL
      
      repunit_trace <- repunit_trace %>% 
        dplyr::mutate(repunit = factor(x = repunit, levels = group_names))  # order repunit
      
    } else {  # groupvec
      
      if(!all(file.exists(paste0(path, "/", mixvec, "_collection_trace.csv")))) {
        
        stop("Not all mixtures in `mixvec` have a `collection_trace.csv` file in `path`!!!")
        
      }  # make sure output files exist
      
      message("    Building trace output from `collection_trace.csv` files, `groupvec`, and `group_names`.")
      
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl, cores = ncores)  
      
      collection_trace <-
        foreach::foreach(mixture = mixvec, .packages = c("tidyverse")) %dopar% {
          
          collection_trace_mix <- suppressMessages(readr::read_csv(file = paste0(path, "/", mixture, "_collection_trace.csv")))
          
          if (!"chain" %in% names(collection_trace_mix)) {
            dplyr::mutate(collection_trace_mix, chain = 1L)
          } # in case older files without multichain
          
          collection_trace_mix <- collection_trace_mix %>% 
            tidyr::pivot_longer(-c(sweep, chain), names_to = "collection", values_to = "pi") %>%
            dplyr::mutate(mixture_collection = mixture)
          
        } %>% dplyr::bind_rows()  # build output from "collection_trace.csv"
      
      parallel::stopCluster(cl)
      
      base_collections <- unique(collection_trace$collection)  # baseline collections are the same order as in rubias output
      
      repunit_new.df <- tibble::tibble(collection = base_collections, repunit_new = group_names[groupvec])  # tibble of new repunit from groupvec
      
      repunit_trace <- collection_trace %>% 
        dplyr::left_join(repunit_new.df, by = "collection") %>%  # join with new repunit
        dplyr::rename(repunit = repunit_new) %>%  # rename new repunit
        dplyr::group_by(mixture_collection, chain, sweep, repunit) %>%  # group and order
        dplyr::summarise(rho = sum(pi), .groups = "drop") # summarise pi (collection) to rho (repunit)
      
    }  # build repunit_trace from "collection_trace.csv", `groupvec`, and `group_names`
    
    # Build `bootstrapped_proportions` if `bias_corr = TRUE`
    if(bias_corr) {
      
      if(!all(file.exists(paste0(path, "/", mixvec, "_bias_corr.csv")))) {
        
        stop("Not all mixtures in `mixvec` have a `bias_corr.csv` file in `path`!!!")
        
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
  ## If rubias_output, ----
  # create `repunit_trace` and `bootstrapped_proportions`
  # only for single chain
  if(!is.null(rubias_output)) {
    
    message("Summarizing results from `rubias output`.")
    
    if(is.null(mixvec)) {
      
      mixvec <- unique(rubias_output$mix_prop_traces$mixture_collection)
      
    }  # used to order mixtures as factor
    
    # Build repunit_trace from `rubias_output` if `groupvec` is NULL
    if(is.null(groupvec)) {
      
      repunit_trace <- rubias_output$mix_prop_traces %>% 
        dplyr::mutate(chain = 1L) %>% 
        dplyr::group_by(mixture_collection, chain, sweep, repunit) %>%  # group to summarize across collections
        dplyr::summarise(rho = sum(pi), .groups = "drop") # summarize collections to repunits
      
      if(is.null(group_names)) {group_names <- unique(repunit_trace$repunit) }  # assign `group_names` from `rubias_output`, if NULL, order may be wrong
      
    } else {
      
      base_collections <- unique(rubias_output$mix_prop_traces$collection)  # baseline collections are the same order as in rubias output
      
      repunit_new.df <- tibble::tibble(collection = base_collections,
                                       repunit_new = group_names[groupvec])  # tibble of new repunit from groupvec
      
      repunit_trace <- rubias_output$mix_prop_traces %>% 
        dplyr::left_join(repunit_new.df, by = "collection") %>%  # join with new repunit
        dplyr::select(-repunit) %>%  # drop old repunit
        dplyr::rename(repunit = repunit_new) %>%  # rename new repunit
        dplyr::group_by(mixture_collection, chain, sweep, repunit) %>%  # group and order
        dplyr::summarise(rho = sum(pi), .groups = "drop")  # summarise pi (collection) to rho (repunit)
      
    }  # build repunit_trace from `rubias_output`, `groupvec`, and `group_names`
    
    # Build `bootstrapped_proportions` if `bias_corr = TRUE`
    if(bias_corr) {
      
      message("    Building bias correction output from `rubias_output`.")
      
      bootstrapped_proportions <- rubias_output$bootstrapped_proportions
      
    }  # bias_corr
    
  }  # rubias_output
  
  #~~~~~~~~~~~~~~~~
  ## Verify that catchvec specifies all mixtures ----
  if(length(mixvec) != length(catchvec)) {
    
    stop("`mixvec` and `catchvec` are not the same length!!!")
    
  }
  
  #~~~~~~~~~~~~~~~~
  ## Calculate `d_rho` for bias correction if specified ----
  if(bias_corr) {
    
    if(nrow(bootstrapped_proportions) == 0) {stop("There is no bias corrected output!!!")}
    
    mixing_proportions_rho <- repunit_trace %>% 
      dplyr::filter(sweep >= burn_in) %>%   # remove burn_in
      dplyr::group_by(mixture_collection, repunit) %>%  # group by mixture and repunit across sweeps
      dplyr::summarise(rho = mean(rho), .groups = "drop") 
    
    d_rho <- mixing_proportions_rho %>% 
      dplyr::left_join(bootstrapped_proportions, by = c("mixture_collection", "repunit")) %>%  # join with `bootstrapped_proportions`
      dplyr::mutate(d_rho = rho - bs_corrected_repunit_ppn) %>%  # calculate d_rho
      dplyr::select(mixture_collection, repunit, d_rho)  # drop other variables
    
  }
  
  #~~~~~~~~~~~~~~~~
  ## Apply bias correction if `d_rho` exists ----
  if(exists("d_rho")) {
    
    repunit_trace <- repunit_trace %>% 
      dplyr::left_join(d_rho, by = c("mixture_collection", "repunit")) %>%  # join trace with d_rho
      dplyr::mutate(rho = rho - d_rho) %>%  # subtract d_rho
      dplyr::select(-d_rho)
    
  }
  
  #~~~~~~~~~~~~~~~~
  ## Roll up to broad-scale groups if `groupvec_new` specified ----
  if(!is.null(groupvec_new)) {
    
    level_key <- sapply(group_names, function(grp) {
      
      i = which(group_names == grp)
      group_names_new[groupvec_new[i]]
      
    }, simplify = FALSE )  # set up level_key to use with recode to roll up groups
    
    repunit_trace <- repunit_trace %>% 
      dplyr::mutate(repunit = dplyr::recode(repunit, !!!level_key)) %>% 
      dplyr::group_by(mixture_collection, chain, sweep, repunit) %>% 
      dplyr::summarise(rho = sum(rho), .groups = "drop") 
    
  }  
  
  #~~~~~~~~~~~~~~~~
  ## Summary statistics for stratified estimate
  
  harvest <- dplyr::tibble(mixture_collection = mixvec, harvest = catchvec)
  
  if(is.null(cv)){
    
    lo_CI = alpha / 2
    hi_CI = 1 - (alpha / 2)
    
    out_sum <- repunit_trace %>% 
      dplyr::filter(mixture_collection %in% mixvec) %>%  # only stratify over mixvec
      dplyr::filter(sweep >= burn_in) %>%  # remove burn_in
      dplyr::left_join(harvest, by = "mixture_collection") %>%  # join harvest data from `catchvec`
      dplyr::mutate(rho_stratified = rho * harvest) %>%  # multiply each strata by strata harvest
      dplyr::group_by(chain, sweep, repunit) %>%  # summarise mixtures by sweep and repunit for each chain
      dplyr::summarise(rho = sum(rho_stratified),
                       rho_pi = sum(rho_stratified)/ sum(harvest),
                       .groups = "keep") %>% # sum up harvest numbers in a stratum by repunit (3/18/22)
      dplyr::mutate(stratified = newname) %>%  # define newname
      dplyr::select(stratified, chain, sweep, repunit, rho, rho_pi) %>% 
      dplyr::group_by(stratified, repunit) %>%  # calculate summary statistics
      dplyr::summarise(mean_harv = mean(rho),
                       sd_harv = sd(rho),
                       median_harv = median(rho),
                       loCI_harv = quantile(rho, probs = lo_CI),
                       hiCI_harv = quantile(rho, probs = hi_CI),
                       mean = mean(rho_pi),
                       sd = sd(rho_pi),
                       median = median(rho_pi),
                       loCI = quantile(rho_pi, probs = lo_CI),
                       hiCI = quantile(rho_pi, probs = hi_CI),
                       `P=0` = mean(rho_pi < threshold), .groups = "drop") %>% # summary statistics to return
      magrittr::set_colnames(c("stratified_mixture", "repunit",
                               "mean_harv", "sd_harv", "median_harv", 
                               paste0(lo_CI * 100, "%_harv"), paste0(hi_CI * 100, "%_harv"),
                               "mean", "sd", "median", 
                               paste0(lo_CI * 100, "%"), paste0(hi_CI * 100, "%"), "P=0"))
  } else {
    
    harvest <- harvest %>% 
      dplyr::mutate(cv = cv)
    
    if(length(mixvec) != length(cv)) {
      stop("`mixvec` and `cv` are not the same length, hoser!!!")
    }
    
    lo_CI <- alpha / 2
    hi_CI <- 1 - (alpha / 2)
    
    if (is.null(group_names_new)) grp_names <- group_names
    else grp_names <- group_names_new
    
    out_sum <- repunit_trace %>% 
      dplyr::filter(mixture_collection %in% mixvec) %>%  # only stratify over mixvec
      dplyr::filter(sweep >= burn_in) %>%  # remove burn_in
      dplyr::full_join(harvest, by = "mixture_collection") %>%
      dplyr::mutate(lnvar = log(cv^2 + 1), lnmean = log(harvest)-log(cv^2 + 1)/2) %>%
      tidyr::pivot_wider(names_from = repunit, values_from = rho) %>%
      dplyr::group_by(mixture_collection) %>%
      dplyr::mutate(h0 = rlnorm(length(unique(sweep))*max(chain), lnmean, sqrt(lnvar))) %>%
      dplyr::ungroup() %>% 
      dplyr::mutate_at(.var = grp_names, .funs = ~.*h0) %>% 
      dplyr::select(-cv, -lnvar, -lnmean) %>% 
      tidyr::pivot_longer(cols = -dplyr::all_of(c("mixture_collection", "chain", "sweep", "harvest", "h0")), names_to = "repunit", values_to = "r_h0") %>% 
      dplyr::mutate(repunit = factor(x = repunit, levels = grp_names)) %>% 
      dplyr::group_by(chain, sweep, repunit) %>% 
      dplyr::summarise(rho = sum(r_h0),
                       rho_pi = sum(r_h0)/sum(h0),
                       tot_harv = sum(h0),
                       .groups = "keep") %>% # sum up harvest numbers in a stratum by repunit (3/18/22)
      dplyr::mutate(stratified = newname) %>%  # define newname
      dplyr::group_by(stratified, repunit) %>% # calculate summary statistics
      dplyr::summarise(mean_harv = mean(rho),
                       sd_harv = sd(rho),
                       median_harv = median(rho),
                       # loCI_harv = HDInterval::hdi(rho, 1- alpha)[1],
                       # hiCI_harv = HDInterval::hdi(rho, 1- alpha)[2],
                       loCI_harv = quantile(rho, probs = lo_CI),
                       hiCI_harv = quantile(rho, probs = hi_CI),
                       mean = mean(rho_pi),
                       sd = sd(rho_pi),
                       median = median(rho_pi),
                       loCI = quantile(rho_pi, probs = lo_CI),
                       hiCI = quantile(rho_pi, probs = hi_CI),
                       `P=0` = mean(rho_pi < threshold),
                       mean_backcalc = mean_harv/mean(tot_harv),
                       sd_backcalc = sd_harv/mean(tot_harv),
                       median_backcalc = median_harv/mean(tot_harv),
                       loCI_backcalc = loCI_harv/mean(tot_harv),
                       hiCI_backcalc = hiCI_harv/mean(tot_harv),
                       .groups = "drop") %>% # summary statistics to return
      dplyr::mutate(hiCI_backcalc = replace(hiCI_backcalc, which(hiCI_backcalc > 1), 1),
                    median_backcalc = replace(median_backcalc, which(median_backcalc > 1), 1),
                    mean_backcalc = replace(mean_backcalc, which(mean_backcalc > 1), 1)) %>% 
      
      magrittr::set_colnames(c("stratified_mixture", "repunit",
                               "mean_harv", "sd_harv", "median_harv", 
                               paste0(lo_CI * 100, "%_harv"), paste0(hi_CI * 100, "%_harv"),
                               "mean", "sd", "median", 
                               paste0(lo_CI * 100, "%"), paste0(hi_CI * 100, "%"), "P=0",
                               "mean_backcalc", "sd_backcalc", "median_backcalc", 
                               paste0(lo_CI * 100, "%_backcalc"),
                               paste0(hi_CI * 100, "%_backcalc")))
    
  }
  
  return(out_sum)
  
}  # end function