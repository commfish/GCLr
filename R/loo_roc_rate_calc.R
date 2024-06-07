#' Calculate calculate assignment rates from leave-one-out analysis
#' 
#' This function takes the leave-one-out output from [rubias::self_assign()] and calculates true positive, false negative, false positive, and true negative assignment error rates for each reporting group. 
#' The output from this function can then be used to create a precision-recall curve plot for determining an appropriate individual assignment threshold and whether a reporting group is sufficiently identifiable for producing individual assignment estimates.
#' 
#' @param data the unmodified output from [rubias::self_assign()]
#' @param thres_levels the assignment thresholds levels that you want to plot
#' @param group_names the reporting group names
#' @param ncores A numeric value for the number of cores to use. 
#' 
#' @seealso [rubias::self_assign()]
#' @seealso [GCLr::plot_loo_prec_rec()]
#' 
#' @note Depending on how large your baseline data set is, this function can take a long time to run.  
#' If you want to be able to work on your computer while the function is running make sure to set `ncores` below the total number of cores on your machine.
#' 
#' @returns This function produces a tibble 10 variables:
#'      \itemize{
#'                \item \code{repunit} (character): reporting group
#'                \item \code{fn} (integer): number of false negative assignments
#'                \item \code{tn} (integer): number of true negative assignments
#'                \item \code{tp} (integer): number of true positive assignments
#'                \item \code{fp} (integer): number of false positive assignments
#'                \item \code{tpr} (double): the true positive assignment rate - tp / (tp + fn)
#'                \item \code{fpt} (double): the false positive assignment rate - fp / (fp + tn)
#'                \item \code{acc} (double): accuracy rate - (tp + tn) / (tp + tn + fp + fn)
#'                \item \code{pre} (double): precision rate - tp / (tp + fp)
#'                \item \code{level} (double): the assignment threshold level (proportion)
#'                \item \code{threshold} (character): the assignment threshold level (percent) for precision-recall curve plots.
#'          }
#' 
#' @examples
#'
#'  Self_Assign <- rubias::self_assign(reference = GCLr::ex_baseline, gen_start_col = 5)
#' 
#'  groups <- GCLr::ex_baseline$repunit %>% unique()
#' 
#'  loo_rate_calc(data = Self_Assign, thres_levels = seq(0.5, .99, by = 0.01), group_names = groups, ncores = 3)
#'  
#' @aliases loo_roc_rate_calc
#' 
#' @rdname loo_rate_calc
#' 
#' @export 
loo_rate_calc <- function(data, thres_levels, group_names, ncores = parallel::detectCores()) {
  
  rubias_loo_vars <- c("indiv", "collection", "repunit", "inferred_collection", "inferred_repunit", 
                       "scaled_likelihood", "log_likelihood", "z_score", "n_non_miss_loci", 
                       "n_miss_loci", "missing_loci")
  
  if(length(setdiff(names(data), rubias_loo_vars)) > 0){
    
    stop("The data object provided does not contain all of the nessecary varibles. 
         Make sure you are using an unmodified object produced by rubias::self_assign().")
     
  }
  
  sa_to_repu <- data %>% 
    dplyr::group_by(indiv, collection, repunit, inferred_repunit) %>%
    dplyr::summarise(repu_scaled_like = sum(scaled_likelihood), .groups = "drop_last")
  
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)
  
  `%dopar%` <- foreach::`%dopar%`
  
  output <- foreach::foreach(repu = group_names, .packages = "tidyverse") %dopar% {
    
    outcomes <- c("tp", "fn", "fp", "tn")
    
      lapply(thres_levels, function(ts) {
        
        counts <- sa_to_repu %>%
          dplyr::filter(inferred_repunit == repu) %>%
          dplyr::mutate(
            outcome = dplyr::case_when(
              inferred_repunit == repunit & repu_scaled_like >= ts ~ "tp", # true positive
              inferred_repunit == repunit & repu_scaled_like < ts ~ "fn", # false negative
              inferred_repunit != repunit & repu_scaled_like >= ts ~ "fp", # false positive
              inferred_repunit != repunit & repu_scaled_like < ts ~ "tn" # true negative
            )
          ) %>%
          dplyr::group_by(outcome) %>%
          dplyr::summarize(n = n(), .groups = "drop") %>% # counts for each outcome
          tidyr::pivot_wider(names_from = outcome, values_from = n)
        
        # Accounting for missing outcomes
        if (sum(outcomes %in% names(counts)) != 4) {
          missing <- outcomes[-match(names(counts), outcomes)]
          
          for (mis in missing) {
            
            counts <- counts %>%
              dplyr::mutate(!!sym(mis) := 0)
            
          }
          
        }
        
        counts %>%
          dplyr::mutate(
            tpr = tp / (tp + fn), # true positive rate (power; sensitivity)
            fpr = fp / (tn + fp), # false positive rate (type 1 error)
            acc = (tp + tn) / (tp + tn + fp + fn), # accuracy rate
            pre = tp / (tp + fp), # precision rate (positive predictive value)
            level = ts
          ) # threshold level
        
      }) %>% bind_rows()
      
    } %>%
    purrr::set_names(group_names) %>% 
    dplyr::bind_rows(.id = "repunit") %>% 
    dplyr::mutate(threshold = paste0(level*100, "%"))
  
  parallel::stopCluster(cl)# Stop cluster
  
  return(output)
  
}
#' @export
loo_roc_rate_calc <- loo_rate_calc