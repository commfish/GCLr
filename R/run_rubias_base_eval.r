#' Run Rubias Baseline Evaluation Tests
#'
#' This function is a wrapper for [GCLr::run_rubias_mix()]. It runs baseline evaluation mixtures produced by [GCLr::create_rubias_base_eval()] and multicores by test group. i.e., mixtures for a test group are run on a different core simultaneously.
#' It is designed for use on a server with many cores and ample RAM.
#'
#' @param tests a tibble produced by [GCLr::base_eval_sample_sizes()] containing `test_group` and `scenario`
#'      
#' @param group_names a character vector of group names, used to sort `repunit` as a factor
#' 
#' @param gen_start_col an integer indicating the starting column index for genetic data in the mixture and baseline files (default: 5)
#' 
#' @param base.path the file path where the baseline .csv files are located
#' 
#' @param mix.path the file path where the mixture .csv files are located
#' 
#' @param out.path the file path to save output from each mixture as a .csv
#' 
#' @param method a choice among "MCMC", "PB", and "BR" methods for estimating mixture proportions (see details)
#' 
#' @param alle_freq_prior a one-element named list specifying the prior to be used when generating Dirichlet parameters for genotype likelihood calculations. Valid methods include "const", "scaled_const", and "empirical". See [rubias::list_diploid_params()] for method details. (default: `list(const_scaled = 1)`)
#' 
#' @param pi_prior The prior to be added to the collection allocations, in order to generate pseudo-count Dirichlet parameters for the simulation of new pi vectors in MCMC. Default value of NA leads to the calculation of a symmetrical prior based on pi_prior_sum. To provide other values to certain collections, you can pass in a data frame with two columns, "collection" listing the relevant collection, and "pi_param" listing the desired prior for that collection (see [GCLr::create_prior()]. Specific priors may be listed for as few as one collection. The special collection name "DEFAULT_PI" is used to set the prior for all collections not explicitly listed; if no "DEFAULT_PI" is given, it is taken to be 1/(# collections).
#' 
#' @param pi_init The initial value to use for the mixing proportion of collections. This lets the user start the chain from a specific value of the mixing proportion vector. If pi_init is NULL (the default) then the mixing proportions are all initialized to be equal. Otherwise, you pass in a data frame with one column named "collection" and the other named "pi_init" (see [GCLr::random_inits()]. Every value in the pi_init column must be strictly positive (> 0), and a value must be given for every collection. If they sum to more than one the values will be normalized to sum to one.
#' 
#' @param reps the number of MCMC iterations (default: 25000)
#' 
#' @param burn_in how many iterations to discard in the beginning of MCMC when doing the mean calculation. They will still be returned in the traces if desired. (default: 5000)
#' 
#' @param pb_iter how many bootstrapped data sets to do for bootstrap correction using method PB. (default: 100)
#' 
#' @param prelim_reps for method "BR", the number of reps of conditional MCMC (as in method "MCMC") to perform prior to MCMC with baseline resampling. The posterior mean of mixing proportions from this conditional MCMC is then used as pi_init in the baseline resampling MCMC (default: NULL)
#' 
#' @param prelim_burn_in for method "BR", this sets the number of sweeps out of prelim_reps that should be discarded as burn in when preparing the posterior means of the mixing proportions to be set as pi_init in the baseline resampling MCMC. (default: NULL)
#' 
#' @param sample_int_Pi how many iterations between storing the mixing proportions trace. Can't be 0. Can't be so large that fewer than 10 samples are taken from the burn in and the sweeps. (default: 10)
#' 
#' @param sample_theta 	for method "BR", whether or not the function should store the posterior mean of the updated allele frequencies. (default: TRUE)
#' 
#' @param pi_prior_sum For pi_prior = NA, the prior on the mixing proportions is set as a Dirichlet vector of length C, with each element being W/C, where W is the pi_prior_sum and C is the number of collections. By default this is 1. If it is made much smaller than 1, things could start to mix more poorly. (default: 1)
#' 
#' @param seed An integer to set the seed for the MCMC (default: 56)
#' 
#' @param ncores A numeric vector of length one indicating the number of cores to use (default: 4)
#' 
#' @param file_type The file type of the baseline and mixture input files, either "fst" or "csv" (default: "fst")(see details).
#'
#' @param out_file_type The file type for the output files, either "fst" or "csv" (default: "fst")(see details).
#'
#' @details "MCMC" estimates mixing proportions and individual posterior probabilities of assignment through Markov-chain Monte Carlo conditional on the reference allele frequencies, while "PB" does the same with a parametric bootstrapping correction, and "BR" runs MCMC sweeps while simulating reference allele frequencies using the genotypes of mixture individuals and allocations from the previous sweep. All methods default to a uniform 1/(# collections or RUs) prior for the mixing proportions.
#'          The function can read in .csv or .fst files.  .fst files are compressed, so they save hard drive space, and they are faster to save and read back into R. .csv is also an option to make the function backwards compatible with older analyzes that produced .csv files.
#' 
#' @return  This function breaks up the \emph{rubias} output by mixture_collection and for each mixture the following are saved as .csv or .fst files:
#'    \itemize{
#'        \item collection level trace, wide format, collections in order of baseline - akin to .BOT file from BAYES
#'        \item repunit level trace, wide format, repunit in order of `group_names` - akin to .RGN file from BAYES
#'        \item straight dump of the `indiv_posteriors` tibble - without column `missing_loci`
#'        \item straight dump of the `indiv_posteriors` tibble - without column `missing_loci`
#'        \item straight dump of the `bootstrapped_proportions`
#'    }
#'    
#' @seealso [GCLr::run_rubias_mix()]
#' @seealso [GCLr::create_rubias_base_eval()]
#'
#' @examples
#' \dontrun{
#' 
#'    sample_sizes <- GCLr::base_eval_sample_sizes(sillyvec = sillyvec, group_names = group_names, groupvec = groupvec, scenarios = round(seq(.01, 1, .01), 2), mixsize = 200, maxprop = 0.5, seed = 123)
#'    GCLr::create_rubias_base_eval(sillyvec = sillyvec, group_names = group_names, test_groups = group_names, loci = loci, groupvec = groupvec, sample_sizes = sample_sizes, prprtnl = TRUE, seed = 123, ncores = parallel::detectCores())
#'    tests <- sample_sizes %>% group_by(test_group, scenario) %>% summarize(test_group = test_group %>% unique(), scenario = scenario %>% unique(), .groups = "drop_last")
#'    GCLr::run_rubias_base_eval(tests = tests, group_names = group_names, gen_start_col = 5, base.path = "rubias/baseline", mix.path = "rubias/mixture", out.path = "rubias/output", seed = 56, ncores = parallel::detectCores())
#' 
#' }
#'    
#' @export
run_rubias_base_eval <- function(tests, group_names, gen_start_col = 5,  base.path = "rubias/baseline", mix.path = "rubias/mixture", out.path = "rubias/output", method = "MCMC", alle_freq_prior = list(const_scaled = 1), pi_prior = NA, 
                                  pi_init = NULL, reps = 25000, burn_in = 5000, pb_iter = 100, prelim_reps = NULL, prelim_burn_in = NULL, sample_int_Pi = 10, sample_theta = TRUE, pi_prior_sum = 1, seed = 56, ncores = parallel::detectCores(), file_type = c("fst", "csv")[1], out_file_type = c("fst", "csv")[1]){
 
  start_time <- Sys.time()
  
  test_groups <- tests$test_group %>% 
    unique()
  
  # Multicore loop
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)
  
  `%dopar%` <- foreach::`%dopar%`
  
  foreach::foreach(g = test_groups, .export = "run_rubias_mix", .packages = c("tidyverse", "rubias", "readr")) %dopar% {
    
    scenario_names <- tests %>% 
      dplyr::filter(test_group==g) %>% 
      tidyr::unite(col = "name", test_group, scenario) %>% 
      dplyr::pull(name)
    
    for(scn in scenario_names){
      
      if(file_type == "fst"){
        
        mixture <- fst::read_fst(path = paste0(mix.path, "/", scn, ".mix.fst")) 
        
        baseline <- fst::read_fst(path = paste0(base.path, "/", scn, ".base.fst"))
        
      }
      
      if(file_type == "csv"){
        
        mixture <- readr::read_csv(file = paste0(mix.path, "/", scn, ".mix.csv"), col_types = cols(.default = "c")) 
        
        baseline <- readr::read_csv(file = paste0(base.path, "/", scn, ".base.csv"), col_types = cols(.default = "c"))
        
      }
      
      GCLr::run_rubias_mix(reference = baseline, 
                         mixture = mixture, 
                         group_names = group_names, 
                         gen_start_col = gen_start_col ,
                         method = method, 
                         alle_freq_prior = alle_freq_prior, 
                         pi_prior = pi_prior, 
                         pi_init = pi_init, 
                         reps = reps, 
                         burn_in = burn_in, 
                         pb_iter = pb_iter, 
                         prelim_reps = prelim_reps, 
                         prelim_burn_in = prelim_burn_in, 
                         sample_int_Pi = 10, 
                         sample_theta = sample_int_Pi, 
                         pi_prior_sum = pi_prior_sum, 
                         file = out.path, 
                         seed = seed,
                         out_file_type = out_file_type)
      
    }
    
  }# end multicore loop
  
  parallel::stopCluster(cl)
  
  message(paste0("Analysis complete!! The test mixture output filese are located in this folder: ", out.path))
  
  Sys.time()-start_time
  
}