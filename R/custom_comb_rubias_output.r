#' Custom Combine \pkg{rubias} Output
#'
#' This function computes summary statistics from \pkg{rubias} output, similar to `CustomCombineBAYESOutput`.
#' However, output is a tibble with `mixture_collection` as a column, instead of each mixture as its own list.
#' It can take either the `rubias_output` list object from `run_rubias_mix()` or [rubias::infer_mixture()],
#' OR it can read in the .csv files created by `run_rubias_mix()`.
#'
#' NOTE: Currently this function only allows bias correction for the reporting groups run in the mixture.
#' It cannot do bias correction for different baseline groupvecs because the current \pkg{rubias} output only
#' gives the bias-corrected means for each `mixture_collection` and `repunit` (i.e., `rho`, not `pi`).
#'
#' UPDATE: This function CAN do bias correction if you are rolling up groups from fine-scale to broad-scale.
#' To use this functionality, specify `group_names` as the original, fine-scale groups, `groupvec_new` as 
#' the groupvec to go from fine-scale to broad-scale groups, and `group_names_new` as the broad-scale groups
#' (i.e., length(groupvec_new) = length(group_names), and max(groupvec_new) == length(group_names_new))
#' Compatible with multichain output.
#'
#' @param rubias_output Output list object from `run_rubias_mix()` or [rubias::infer_mixture()].
#' @param mixvec Character vector of mixture sillys, used to read in output .csv files if `rubias_output = NULL`.
#' @param group_names Character vector of group_names, used to sort repunit as a factor, can get from .csv.
#' @param group_names_new Character vector of new group_names, used to roll up groups from fine-scale to broad-scale for bias correction
#'                       if specified, `groupvec_new` must be = length(group_names_new).
#' @param groupvec Numeric vector indicating the group affiliation of each pop in sillyvec, used if resummarizing to new groups.
#' @param groupvec_new A numeric vector indicating the new group affiliation of each group, used if resummarizing fine-scale groups to broad-scale groups with bias correction.
#' @param path Character vector of where to find output from each mixture as a .csv (created by `run_rubias_mix()`).
#' @param alpha Numeric constant specifying credibility intervals, default is 0.1, which gives 90% CIs (i.e., 5% and 95%).
#' @param burn_in Numeric constant specifying how many sweeps were used for burn_in in `run_rubias_mix()` or [rubias::infer_mixture()].
#' @param bias_corr Logical switch indicating whether you want bias-corrected values from `method = "PB"` or not, 
#'                 currently can NOT do bias correction if not using the same repunits that were run in the mixture.
#' @param threshold Numeric constant specifying how low stock comp is before assuming 0, used for `P=0` calculation, default is from BAYES.
#' @param plot_trace Logical switch, when on will create a trace plot for each mixture and repunit (reporting group).
#' @param ncores A numeric vector of length one indicating the number of cores to use (ncores is only used when is.null(rubias_output) == TRUE).
#'
#' @return A tibble with 8 fields for each mixture and repunit (reporting group).
#'   - \code{mixture_collection}: Factor of mixtures (only a factor for ordering, plotting purposes).
#'   - \code{repunit}: Factor of reporting groups (only a factor for ordering, plotting purposes).
#'   - \code{mean}: Mean stock composition.
#'   - \code{sd}: Standard deviation.
#'   - \code{median}: Median stock composition.
#'   - \code{loCI}: Lower bound of credibility interval.
#'   - \code{hiCI}: Upper bound of credibility interval.
#'   - \code{P=0}: The proportion of the stock comp distribution that was below `threshold` (i.e., posterior probability that stock comp = 0).
#'
#' @seealso custom_comb_bayes_output()
#'
#' @examples
#' \dontrun{
#' load("V:/Analysis/1_SEAK/Sockeye/Mixture/Lynn Canal Inseason/2018/OLD rubias/output/test/custom_comb_rubias_output_test.RData")
#' lynncanal_2015_SW26_27.sum <- custom_comb_rubias_output(rubias_output = lynncanal_test.out, group_names = LynnCanal_groups7, bias_corr = TRUE)
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

  # Error catching ----
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

  if(!is.null(groupvec_new)) {
    if(is.null(group_names)) {
      stop("Need to provide `group_names` if introducing a new `groupvec`, hoser!!!")
    }
    if(length(groupvec_new) != length(group_names) | max(groupvec_new) != length(group_names_new))
      stop("If specifying `groupvec_new`, you must be rolling up from fine-scale groups to broad-scale groups (i.e. length(groupvec_new) == length(group_names), and max(groupvec_new) == length(group_names_new))")
  }
  
  `%dopar%` <- foreach::`%dopar%`
  
  # Summaries ----
  #~~~~~~~~~~~~~~~~
  ## If no rubias_output, ----
  # make from .csv files
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
      
    } else {  # groupvec
      
      if(!all(file.exists(paste0(path, "/", mixvec, "_collection_trace.csv")))) {
        stop("Not all mixtures in `mixvec` have a `collection_trace.csv` file in `path`, hoser!!!")
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
      repunit_new.df <- tibble::tibble(collection = base_collections,
                                       repunit_new = group_names[groupvec])  # tibble of new repunit from groupvec
      repunit_trace <- collection_trace %>% 
        dplyr::left_join(repunit_new.df, by = "collection") %>%  # join with new repunit
        dplyr::rename(repunit = repunit_new) %>%  # rename new repunit
        dplyr::group_by(mixture_collection, chain, sweep, repunit) %>%  # group and order
        dplyr::summarise(rho = sum(pi), .groups = "drop")   # summarise pi (collection) to rho (repunit)
        
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
      
      repunit_trace <-
        rubias_output$mix_prop_traces %>%
        dplyr::mutate(
          chain = if ("chain" %in% names(rubias_output$mix_prop_traces)) {
            chain
            } else {1L}
          ) %>% 
        dplyr::group_by(mixture_collection, chain, sweep, repunit) %>%  # group to summarize across collections
        dplyr::summarise(rho = sum(pi), .groups = "drop" ) # summarize collections to repunits
        
      if(is.null(group_names)) {group_names <- unique(repunit_trace$repunit)}  # assign `group_names` from `rubias_output`, if NULL, order may be wrong
      
    } else {
      
      base_collections <- unique(rubias_output$mix_prop_traces$collection)  # baseline collections are the same order as in rubias output
      repunit_new.df <- tibble::tibble(collection = base_collections,
                                       repunit_new = group_names[groupvec])  # tibble of new repunit from groupvec
      repunit_trace <-
        rubias_output$mix_prop_traces %>% 
        dplyr::mutate(
          chain = if ("chain" %in% names(rubias_output$mix_prop_traces)) {
            chain
          } else {1L}
        ) %>% 
        dplyr::left_join(repunit_new.df, by = "collection") %>%  # join with new repunit
        dplyr::select(-repunit) %>%  # drop old repunit
        dplyr::rename(repunit = repunit_new) %>%  # rename new repunit
        dplyr::group_by(mixture_collection, chain, sweep, repunit) %>%  # group and order
        dplyr::summarise(rho = sum(pi), .groups = "drop")   # summarise pi (collection) to rho (repunit)
        
    }  # build repunit_trace from `rubias_output`, `groupvec`, and `group_names`
    
    # Build `bootstrapped_proportions` if `bias_corr = TRUE`
    if(bias_corr) {
      message("    Building bias correction output from `rubias_output`.")
      bootstrapped_proportions <- rubias_output$bootstrapped_proportions
    }  # bias_corr
    
  }  # rubias_output
  
  #~~~~~~~~~~~~~~~~
  ## Calculate `d_rho` for bias correction if specified ----
  if(bias_corr) {
    if(nrow(bootstrapped_proportions) == 0) {stop("There is no bias corrected output, hoser!!!")}
    
    mixing_proportions_rho <- repunit_trace %>% 
      dplyr::filter(sweep >= burn_in) %>%   # remove burn_in
      dplyr::group_by(mixture_collection, repunit) %>%  # group by mixture and repunit across sweeps
      dplyr::summarise(rho = mean(rho), .groups = "drop") 
    
    d_rho <- mixing_proportions_rho %>% 
      dplyr::left_join(bootstrapped_proportions, by = c("mixture_collection", "repunit")) %>%  # join with `bootstrapped_proportions`
      dplyr::mutate(d_rho = rho - bs_corrected_repunit_ppn) %>%  # calculate d_rho
      dplyr::select(mixture_collection, repunit, d_rho)  # drop other variables
  }
  
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
    grp_names <- group_names_new  # for factoring repunit
  } else {
    grp_names <- group_names  # for factoring repunit
  } 
  
  #~~~~~~~~~~~~~~~~
  ## Plot repunit trace ----
  if (plot_trace == TRUE) {
    if (max(repunit_trace$chain) == 1) {
      trace_plot <- repunit_trace %>% 
        dplyr::mutate(
          mixture_collection = factor(x = mixture_collection, levels = mixvec),
          repunit = factor(x = repunit, levels = grp_names)) %>%  # order mixture_collection, repunit
        ggplot2::ggplot(ggplot2::aes(x = sweep, y = rho, colour = repunit)) +
        ggplot2::geom_line() +
        ggplot2::ylim(0, 1) +
        ggplot2::geom_vline(xintercept = burn_in) +
        # ggplot2::annotate(geom = "text", x = burn_in / 2, y = 0.9, label = "Burn-in") +
        ggplot2::theme(legend.position = "none",
                       axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggplot2::facet_grid(repunit ~ mixture_collection)
      
    } else if (max(repunit_trace$chain) > 1) {
      trace_plot <- repunit_trace %>% 
        dplyr::mutate(
          mixture_collection = factor(x = mixture_collection, levels = mixvec),
          repunit = factor(x = repunit, levels = grp_names),
          chain = as.factor(chain)
          ) %>%  # order mixture_collection, repunit
        ggplot2::ggplot(ggplot2::aes(x = sweep, y = rho, colour = chain)) +
        ggplot2::geom_line() +
        # ggplot2::ylim(0, 1) +
        ggplot2::geom_vline(xintercept = burn_in) +
        # ggplot2::annotate(geom = "text", x = burn_in / 2, y = 0.9, label = "Burn-in") +
        ggplot2::theme(legend.position = "none",
                       axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggplot2::facet_grid(repunit ~ mixture_collection)
    } # multi chain
  } # trace plot

  #~~~~~~~~~~~~~~~~
  ## Summary statistics ----
  loCI = alpha / 2
  hiCI = 1 - (alpha / 2)
  nchains <- max(repunit_trace$chain)
  
  mcmc_tr_mc <-
    lapply(unique(repunit_trace$mixture_collection), function(mix) {
      lapply(unique(repunit_trace$repunit), function(grp) {
        sapply(unique(repunit_trace$chain), function(ch) {
          dplyr::filter(repunit_trace,
                        mixture_collection == mix, repunit == grp,
                        chain == ch, sweep >= burn_in) %>%
            dplyr::select(rho) %>%
            dplyr::mutate(rho = coda::mcmc(rho))
        }) %>% coda::as.mcmc.list(.)
      })
    })
  
  out_sum <- repunit_trace %>% 
    dplyr::filter(sweep >= burn_in) %>%  # remove burn_in
    dplyr::mutate(
      mixture_collection = factor(x = mixture_collection, levels = mixvec),
      repunit = factor(x = repunit, levels = grp_names)
      ) %>%  # order mixture_collection, repunit
    dplyr::group_by(mixture_collection, repunit) %>%
    dplyr::summarise(mean = mean(rho),
                     sd = sd(rho),
                     median = median(rho),
                     loCI = quantile(rho, probs = loCI),
                     hiCI = quantile(rho, probs = hiCI),
                     `P=0` = sum(rho < threshold) / length(rho),
                     .groups = "drop") %>%    # summary statistics to return
    dplyr::mutate(
      GR = {if (nchains > 1) {
        lapply(seq(length(mixvec)), function(m) {
          lapply(seq(length(grp_names)), function(g) {
            coda::gelman.diag(mcmc_tr_mc[[m]][[g]],
                              transform = TRUE,
                              autoburnin = FALSE,
                              multivariate = FALSE)$psrf[,"Point est."]
          }) %>% dplyr::bind_rows()
        }) %>% dplyr::bind_rows() %>% unlist()
      } else {NA}},
      n_eff =
        lapply(seq(length(mixvec)), function(m) {
          lapply(seq(length(grp_names)), function(g) {
            coda::effectiveSize(mcmc_tr_mc[[m]][[g]])
            }) %>% dplyr::bind_rows()
        }) %>% dplyr::bind_rows() %>% unlist()
    ) %>%
    dplyr::mutate(loCI = replace(loCI, which(loCI < 0), 0),
                  hiCI = replace(hiCI, which(hiCI < 0), 0),
                  median = replace(median, which(median < 0), 0),
                  mean = replace(mean, which(mean < 0), 0)) %>% 
    dplyr::mutate(loCI = replace(loCI, which(loCI > 1), 1),
                  hiCI = replace(hiCI, which(hiCI > 1), 1),
                  median = replace(median, which(median > 1), 1),
                  mean = replace(mean, which(mean > 1), 1)) %>% 
    magrittr::set_colnames(c("mixture_collection", "repunit", "mean", "sd", "median", 
                             paste0(loCI * 100, "%"), paste0(hiCI * 100, "%"), "P=0",
                             "GR", "n_eff"))
  
  if(exists("trace_plot")) {suppressWarnings(print(trace_plot))}  # plot the trace for each mixture, repunit
  
  return(out_sum)
  }

