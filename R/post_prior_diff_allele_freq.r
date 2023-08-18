#' Get Difference in Allele Frequencies Between Baseline and BAYES Posteriors
#' 
#' This function reads in BAYES MCMC samples of baseline frequencies (.FRQ output files) for each population (one file per chain) and compares them to the original baseline allele frequencies (see details).
#' 
#' @param sillyvec a character vector of population sillys used to create the baseline file.
#' 
#' @param group_names a character vector of reporting group names, where \code{length(group_names) == max(`groupvec`)}
#' 
#' @param groupvec a numeric vector indicating the reporting group affiliation of each population in `sillyvec`, where \code{length(groupvec) == length(sillyvec)}.
#' 
#' @param loci Character vector of the loci used to produce the baseline and mixture files.
#' 
#' @param dir the director where the BAYES mixture output files are located.  e.g., "BAYES/output" This directory should contain a folder named `mixname`, where the .FRQ files are saved.
#' 
#' @param nchains the number of MCMC chains analyzed for the mixture
#' 
#' @param mixname the name of the mixture. There should be a folder in `dir` named the same as the mixture
#' 
#' @param burn The proportion of iterations to drop from the beginning of each chain.
#'                For example, for 40,000 iterations setting burn = 0.5 will drop the first 20,000 iterations.
#' 
#' @param ncores the number of cores to use when [GCLr::calc_freq_pop()] is called on (default = 4)
#' 
#' @details 
#' During a mixed stock analysis in BAYES, the original allele frequencies get updated at each MCMC iteration. 
#' Looking at the difference between the average MCMC allele frequency and the original baseline frequencies can be useful when looking into reasons for non-convergence among chains (i.e., Gelman-Ruban shrink factor > 1.2) 
#' Populations/reporting groups with large changes in allele frequencies for certain population/reporting groups can indicate that the mixture contains fish from a distinct population that is not represented in the baseline.
#'
#' @note To run this function, you must have the baseline silly (.gcl objects) and LocusControl loaded in your current workspace and have BAYES frequency files for each chain of the mixture you are interested in.
#' 
#' @seealso See [BAYES manual](system.file("BAYES", "MANUAL.DOC", package = "GCLr")) for additional details.
#' 
#' @returns a list containing 3 elements:
#'   \itemize{
#'     \item \code{MeanAbsDiffByPopLocus}: a tibble with 4 variables: 
#'     \describe{
#'       \item{\code{chain}}{character; the chain number (e.g. Chain1, Chain2, Chain3, etc.)}
#'       \item{\code{pop}}{factor; the population silly}
#'       \item{\code{locus}}{factor; the locus name}
#'       \item{\code{mean_abs_diff}}{double; the mean absolute difference between the original baseline frequency and the MCMC allele frequencies at each iteration}
#'       }
#'    \item \code{MeanAbsDiffByPop}: a tibble with 4 variables: 
#'     \describe{
#'       \item{\code{chain}}{character; the chain number (e.g. Chain1, Chain2, Chain3, etc.)}
#'       \item{\code{pop}}{factor; the population silly}
#'       \item{\code{mean_abs_diff}}{double; the mean absolute difference between the original population baseline frequency and the group MCMC allele frequencies at each iteration over all loci}
#'       \item{\code{outlier}}{logical; indicates whether the `mean_abs_diff` is an outlier}
#'       }
#'    \item \code{MeanAbsDiffByGroup}: a tibble with 4 variables: 
#'     \describe{
#'       \item{\code{chain}}{character; the chain number (e.g. Chain1, Chain2, Chain3, etc.)}
#'       \item{\code{group}}{character; the group name}
#'       \item{\code{mean_abs_diff}}{double; the mean absolute difference between the original group baseline frequency and the group MCMC allele frequencies at each iteration over all loci}
#'       \item{\code{outlier}}{logical; indicates whether the `mean_abs_diff` is an outlier}
#'       }
#' }
#'   
#' @examples
#' \dontrun{
#' dir <- "V:/Analysis/2_Central/Sockeye/Cook Inlet/Missing Baseline Analysis/BAYES/output"
#' load("V:/Analysis/2_Central/Sockeye/Cook Inlet/Missing Baseline Analysis/Missing Baseline Analysis.RData")
#' post_prior_diff_allele_freq(sillyvec = PopNames69, group_names = Groups, groupvec = groupvec69, loci = loci, dir = dir, nchains = 7, mixname = "Western06early")
#' 
#' }
#' 
#' @export
post_prior_diff_allele_freq <- function(sillyvec, group_names, groupvec, loci, dir, nchains, mixname, burn = 0.5, ncores = 4){
  
  baseline <- GCLr::calc_freq_pop(sillyvec, loci, ncores)
  
  nalleles <- baseline %>% 
    dplyr::group_by(locus) %>% 
    dplyr::summarize(alleles = max(allele_no)) %>% 
    tibble::deframe()
  
  C <- length(sillyvec)
  
  L <- length(loci)
  
  ChainNames <- paste0("Chain", 1:nchains)
  
  mydirs <- paste0(dir, "/", mixname, "/", mixname, ChainNames, "FRQ.FRQ")
  
  myfiles <- lapply(mydirs, function(dr){
    
    myfile <- readr::read_lines(dr)
    
    myfile[(length(myfile)*burn+1):length(myfile)] %>%
      tibble::as_tibble() %>% 
      tidyr::separate(col = value, into = c(NA, "iteration","pop", "locus", "n", paste0("allele", 1:max(nalleles))), sep = "\\s{1,}")
    
  }) %>% 
    setNames(ChainNames) %>% 
    dplyr::bind_rows(.id = "chain")
  
  myQ <- myfiles %>% 
    dplyr::mutate(pop = factor(sillyvec[as.numeric(pop)], levels = sillyvec),
           locus = factor(loci[as.numeric(locus)], levels = loci)) %>% 
    tidyr::pivot_longer(dplyr::all_of(paste0("allele", 1:max(nalleles))), names_to = "allele_no", values_to = "updated_proportion") %>% 
    dplyr::mutate(allele_no = gsub(pattern = "allele", replacement = "", x = allele_no) %>% as.numeric(),
           updated_proportion = as.numeric(updated_proportion)) %>% 
    dplyr::group_by(chain, pop, locus, allele_no) %>% 
    dplyr::summarize(mean_updated_proportion = mean(updated_proportion), .groups = "drop") %>% 
    dplyr::left_join(baseline, by = c("pop" = "silly", "locus", "allele_no")) 
  
  MeanAbsDiffByPopLocus <- myQ %>% 
    dplyr::filter(allele_no == 1) %>% 
    dplyr::mutate(pop = factor(pop, levels = sillyvec),
           locus = factor(locus, levels = loci)) %>% 
    dplyr::group_by(chain, pop, locus) %>% 
    dplyr::summarize(mean_abs_diff = abs(mean_updated_proportion-proportion), .groups = "drop") 
  
  MeanAbsDiffByPop <- MeanAbsDiffByPopLocus %>% 
    dplyr::group_by(chain, pop) %>% 
    dplyr::summarize(mean_abs_diff = mean(mean_abs_diff), .groups = "drop") %>% 
    dplyr::group_by(chain) %>% 
    dplyr::mutate(outlier = outliers::outlier(mean_abs_diff, logical = TRUE))
  
  MeanAbsDiffByGroup <- MeanAbsDiffByPop %>% 
    dplyr::left_join(tibble::tibble(pop = sillyvec, group = group_names[groupvec]), by = "pop") %>% 
    dplyr::group_by(chain, group) %>% 
    dplyr::summarize(mean_abs_diff = mean(mean_abs_diff), .groups = "drop") %>% 
    dplyr::group_by(chain) %>% 
    dplyr::mutate(outlier = outliers::outlier(mean_abs_diff, logical = TRUE))
  
  return(list(MeanAbsDiffByPopLocus = MeanAbsDiffByPopLocus, MeanAbsDiffByPop = MeanAbsDiffByPop, MeanAbsDiffByGroup = MeanAbsDiffByGroup))
  
}