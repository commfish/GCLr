#' Run rubias Mixture Analysis
#'
#' This function is a wrapper for [rubias::infer_mixture()]. It performs the \pkg{rubias} mixture analysis using the provided reference and mixture data frames. 
#' The function estimates population proportions and allele frequencies for the mixture samples given the reference samples.
#'
#' @param reference A data frame containing the reference samples (i.e., baseline object).
#' @param mixture A data frame containing the mixture samples.
#' @param group_names A character vector specifying the names of reference groups.
#' @param gen_start_col An integer specifying the starting column for genotype data in the data frames.
#' @param method The method used for estimation. Default is "MCMC".
#' @param alle_freq_prior A list with a single element named 'const_scaled', specifying the constant scaling for the prior on allele frequencies. Default is 1.
#' @param pi_prior A data frame or tibble with two variables 'collection' and 'pi_param', specifying the prior information for mixture proportions. Default is NA.
#' @param pi_init A matrix or data frame with initial values for mixture proportions. Default is NULL.
#' @param reps An integer specifying the number of MCMC repetitions. Default is 25000.
#' @param burn_in An integer specifying the number of MCMC burn-in iterations. Default is 5000.
#' @param pb_iter An integer specifying the number of parametric bootstrap iterations. Default is 100.
#' @param prelim_reps An integer specifying the number of preliminary repetitions for estimating mixture proportions. Default is NULL.
#' @param prelim_burn_in An integer specifying the number of burn-in iterations for preliminary runs. Default is NULL.
#' @param sample_int_Pi An integer specifying the interval to save MCMC samples of mixture proportions. Default is 10.
#' @param sample_theta A logical value indicating whether to save MCMC samples for allele frequencies. Default is TRUE.
#' @param pi_prior_sum A numeric value specifying the sum constraint for the Dirichlet prior on mixture proportions. Default is 1.
#' @param file A character string representing the file path to save output. Default is "rubias/output".
#' @param seed An integer specifying the random seed for reproducibility. Default is 56.
#' @param nchains Run multiple chains for \pkg{rubias}. Default = 1 for running single `rubias`.
#'
#' @return A list containing the results of the \pkg{rubias} mixture analysis, including estimated mixture proportions, allele frequencies, and other relevant information.
#'
#' @details This function performs the \pkg{rubias} mixture analysis by estimating the proportions of each reference group in the mixture samples 
#' and the allele frequencies of each reference group. It uses the provided method for estimation, such as Markov Chain Monte Carlo (MCMC) 
#' or parametric bootstrap (PB) for uncertainty quantification. The output is saved as .csv files for further analysis and visualization.
#'
#' The `reference` and `mixture` data frames should have the same structure with matching variable names. The `group_names` should represent the 
#' names of reference groups, and these names should correspond to the unique values in the `repunit` column of the `reference` data frame.
#'
#' The `pi_prior` argument allows users to specify prior information for mixture proportions. It should be a data frame or tibble with two variables, 
#' 'collection' and 'pi_param', representing the collection name and prior parameter for the mixture proportions, respectively.
#'
#' When the method is set to "PB" (parametric bootstrap), the function also saves parametric bootstrap bias corrections in separate .csv files.
#'
#' @seealso \code{\link{rubias::infer_mixture}}
#'
#' @examples
#'
#' dir.create(path = path.expand("~/rubias/output"), recursive = TRUE)
#' 
#' group_names <- GCLr::ex_baseline$repunit %>% unique()
#' 
#' run_rubias_mix(reference = GCLr::ex_baseline, mixture = GCLr::ex_mixtures, group_names = group_names, gen_start_col = 5, file = path.expand("~/rubias/output"))
#' 
#' @export
run_rubias_mix <- function(reference, mixture, group_names, gen_start_col, method = "MCMC",
                           alle_freq_prior = list(const_scaled = 1), pi_prior = NA, 
                           pi_init = NULL, reps = 25000, burn_in = 5000, pb_iter = 100,
                           prelim_reps = NULL, prelim_burn_in = NULL,
                           sample_int_Pi = 10, sample_theta = TRUE, pi_prior_sum = 1,
                           file = "rubias/output", seed = 56, nchains = 1) {
  
  if(!dir.exists(file)) {stop("the file path to save output does not exist!")}
  
  if(sum(names(reference)!=names(mixture))>0){
    
    stop("The reference and mixture data frames differ in structure; check number of columns and variable names. 
         Are you using an old rubias reference object with a new mixture object? 
         Old reference objects may have locus headers with periods replacing hyphens.")
    
  }
  
  if(!(all(group_names %in% unique(reference$repunit)) & all(unique(reference$repunit) %in% group_names))) {
    
    stop("Mismatch between `group_names` and `reference`")
    
  }
  
  if(!is.na(pi_prior) %>% as.vector() %>% unique()){
    
    if(!is.data.frame(pi_prior)|sum(names(pi_prior) %in% c("collection", "pi_param")) < 2 ){
      
      stop("pi_prior must be data frame or tibble with two variables 'collection' and 'pi_param'")
      
    }
    
  }
  
  if (method == "PB" & nchains > 1) {
    nchains <- 1L
  }
  
  # Run infer mixture ----
  if (nchains == 1) {
    set.seed(seed = seed)
    
    rubias_out <-
      rubias::infer_mixture(
        reference = reference,
        mixture = mixture,
        gen_start_col = gen_start_col,
        method = method,
        alle_freq_prior = alle_freq_prior,
        pi_prior = pi_prior,
        pi_init = NULL,
        reps = reps,
        burn_in = burn_in,
        pb_iter = pb_iter,
        prelim_reps = prelim_reps,
        prelim_burn_in = prelim_burn_in,
        sample_int_Pi = sample_int_Pi,
        sample_theta = sample_theta,
        pi_prior_sum = pi_prior_sum
      )
    
    rubias_out$mix_prop_traces <-
      rubias_out$mix_prop_traces %>%
      dplyr::mutate(chain = 1)
    
    rubias_out$indiv_posteriors <-
      rubias_out$indiv_posteriors %>% 
      dplyr::select(-missing_loci)  # remove this unnecessary list object
    
  } else {
    chains <- seq(nchains)
    cl <- parallel::makePSOCKcluster(nchains)
    doParallel::registerDoParallel(cl, cores = nchains)
    doRNG::registerDoRNG(seed, once = TRUE)
    `%dorng%` <- doRNG::`%dorng%`
    
    rubias_out00 <- foreach::foreach(
      ch = chains, .packages = c("rubias")
    ) %dorng% {
      rubias::infer_mixture(
        reference = reference,
        mixture = mixture,
        gen_start_col = gen_start_col,
        method = method,
        alle_freq_prior = alle_freq_prior,
        pi_prior = pi_prior,
        pi_init = pi_init[[ch]],
        reps = reps,
        burn_in = burn_in,
        # pb_iter = pb_iter,
        prelim_reps = prelim_reps,
        prelim_burn_in = prelim_burn_in,
        sample_int_Pi = sample_int_Pi,
        sample_theta = sample_theta,
        pi_prior_sum = pi_prior_sum)
      } # dorng
    
    parallel::stopCluster(cl)
    
    rubias_out <- list()
    
    rubias_out$mix_prop_traces <-
      lapply(1:nchains, function(i) {
        dplyr::mutate(rubias_out00[[i]]$mix_prop_traces, chain = i)
      }) %>%
      dplyr::bind_rows()
    
    # indiv_posteriors are the mean of posterior means
    # (prob not the best way to do this)
    # the rest of values are the same for all chains

    rubias_out$indiv_posteriors <- 
      lapply(1:nchains, function(i) {
        rubias_out00[[i]]$indiv_posteriors
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(mixture_collection, indiv, repunit, collection,
                      log_likelihood, z_score, n_non_miss_loci, n_miss_loci) %>%
      dplyr::summarise(PofZ = mean(PofZ), .groups = "drop")
    
  } # else
  
  # Save output ----
  message("Saving output as .csv files")
  mix_sillys <- unique(mixture$collection)
  baseline_pops <- unique(reference$collection)  # correctly ordered
  
  ## Save mix_prop_traces ----
  ### Save at collection level ----
  message("  saving collection traces.", appendLF = FALSE)
  
  time_coll_trace <- system.time({
    invisible(sapply(mix_sillys, function(mixture) {
      
      mix_prop_trace_wide_pi <- rubias_out$mix_prop_traces %>%
        dplyr::filter(mixture_collection == mixture) %>%  # filter to mixture
        dplyr::mutate(collection = factor(x = collection, levels = baseline_pops)) %>%  # use factor to order collections same as baseline
        dplyr::select(sweep, chain, collection, pi) %>%
        tidyr::pivot_wider(names_from = collection, values_from = pi)  # make wide
      
      readr::write_csv(x = mix_prop_trace_wide_pi, file = paste0(file, "/", mixture, "_collection_trace.csv"))
      
    }))
    
  })
  message("   time: ", sprintf("%.2f", time_coll_trace["elapsed"]), 
          " seconds")
  
  ### Save at repuinit level ----
  message("  saving repunit traces.", appendLF = FALSE)
  
  time_repunit_trace <- system.time({
    invisible(sapply(mix_sillys, function(mixture) {
      
      mix_prop_trace_wide_rho <- rubias_out$mix_prop_traces %>%
        dplyr::filter(mixture_collection == mixture) %>%  # filter to mixture
        dplyr::mutate(repunit = factor(x = repunit, levels = group_names)) %>%  # use factor to order repunit same as group_names
        dplyr::group_by(sweep, chain, repunit) %>% 
        dplyr:: summarise(rho = sum(pi), .groups = "drop") %>% 
        dplyr::select(sweep, chain, repunit, rho) %>%
        tidyr::pivot_wider(names_from = repunit, values_from = rho)  # make wide
      
      readr::write_csv(x = mix_prop_trace_wide_rho, file = paste0(file, "/", mixture, "_repunit_trace.csv"))
      
    }))
    
  })
  
  message("   time: ", sprintf("%.2f", time_repunit_trace["elapsed"]), 
          " seconds")
  
  ## Save indiv_posteriors ----
  message("  saving individual posteriors.", appendLF = FALSE)
  
  time_indiv_posteriors <- system.time({
    invisible(sapply(mix_sillys, function(mixture){
      
      indiv_posteriors <- rubias_out$indiv_posteriors %>%
        dplyr::filter(mixture_collection == mixture)
      
      readr::write_csv(x = indiv_posteriors, file = paste0(file, "/", mixture, "_indiv_posteriors.csv"))
    }))
    
  })
  
  message("   time: ", sprintf("%.2f", time_indiv_posteriors["elapsed"]), 
          " seconds")
  
  ## Save bootstrapped_proportions ----
  if(method == "PB") {
    
    message("  saving parametric bootstrap bias corrections.", appendLF = FALSE)
    
    time_bias <- system.time({
      invisible(sapply(mix_sillys, function(mixture) {
        
        bias_corr <- rubias_out$bootstrapped_proportions %>%
          dplyr::filter(mixture_collection == mixture)
        
        readr::write_csv(x = bias_corr, file = paste0(file, "/", mixture, "_bias_corr.csv"))
        
      } ))
    })
    message("   time: ", sprintf("%.2f", time_bias["elapsed"]), 
            " seconds")
    message("Run method = 'PB' with defualt nchains = 1. It is the way.")
  }  # PB
  
  return(rubias_out)
  
}