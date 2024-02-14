#' Get assignment threshold by reporting group
#' 
#' This function takes the unmodified output from `rubias::self_assign()`, 
#' calculates true positive rate and precision, and returns a tibble of minimum
#' and maximum assignment thresholds for filtering top individual assignment results 
#' in order to meet the desired true positive rate and precision.
#' 
#' @param input the unmodified output object from `rubias::self_assign`
#' 
#' @param tpr_cutoff the true positive rate cutoff; default = 0.80
#' 
#' @param pre_cutoff the precision cutoff; default = 0.95
#' 
#' @param min_thres minimum assignment threshold to include in the output; default = 0.50
#' 
#' @param thres_levels a vector of assignment thresholds; `default = seq(0.01, 0.99, by = 0.01)`
#' 
#' @details
#' In an individual assignment analysis, each mixture samples is given a likelihood of assigning to each baseline reporting group. 
#' The assignments for each sampled are filtered for the reporting group with the highest assignment (top assignment). An assignment threshold 
#' is then used to filter out samples with lower top assignments to reduce the number of incorrect assignments. This function takes the leave-one-out 
#' results from `rubias::self_assign()`and summarizes the number of true positives (TP), false negatives (FN), false positives (FP) at different assignment thresholds,
#' and calculates the true positive rate (TPR) and Precision at each threshold and returns a tibble of minimum and maximum assignment thresholds for
#' each reporting group where `TPR >= tpr_cutoff`, `Precision >= pre_cutoff`, and the `assignment threshold >= min_thres`.
#' 
#' \itemize{ 
#'          \item \eqn{TPR = TP / (TP + FN)}
#'          \item \eqn{Precision = TP / (TP + FP)}
#'          }
#' 
#' @returns a tibble with the re
#' 
#' @examples
#' input <- rubias::self_assign(reference = GCLr::ex_baseline, gen_start_col = 5)
#' 
#' assignment_thresholds(input = input, tpr = 0.80, pre = 0.95, min_thres = 0.50, thres_levels = seq(0.01, 0.99, by = 0.01))
#' 
#' @export
IA_thresholds <- function(input, tpr_cutoff = 0.80, pre_cutoff = 0.95, min_thres = 0.50, thres_levels = seq(0.01, 0.99, by = 0.01)){
  
  group_names <- input$repunit %>% 
    unique()
  
  rates <- GCLr::loo_roc_rate_calc(data = input, thres_levels = seq(0.01, 0.99, by = 0.01), group_names = group_names, ncores = parallel::detectCores())
    
  output <- rates %>% 
      dplyr::filter(level >= min_thres, tpr >= tpr_cutoff, pre >= pre_cutoff) %>% 
      dplyr::group_by(repunit) %>% 
      dplyr::summarize(min_threshold = min(level), max_threshold = max(level))
  
  return(output)

  }