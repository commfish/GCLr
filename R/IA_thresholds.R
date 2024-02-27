#' Get assignment threshold by reporting group
#' 
#' This function takes the unmodified output from [rubias::self_assign], 
#' calculates recall (aka true positive rate) and precision, and returns a tibble of minimum
#' and maximum assignment thresholds for filtering top individual assignment results 
#' in order to meet the desired true positive rate and precision.
#' 
#' @param sa_input the unmodified output object from  [rubias::self_assign]
#' 
#' @param rates optional; the unmodified output object from [GCLr::loo_roc_rate_calc] created using the supplied `sa_input`. If left `NULL`, the rates will be calculated from `sa_input`. Supplying a rates object will make this function run a lot faster.
#'
#' @param min_rec the recall cutoff; default = 0.80
#' 
#' @param min_pre the precision cutoff; default = 0.95
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
#' and calculates the recall and precision at each threshold and returns a tibble of minimum and maximum assignment thresholds for
#' each reporting group where `Recall >= min_rec`, `Precision >= min_pre`, and the `assignment threshold >= min_thres`. 
#' TPR and precision are are calculated from the loo results with [GCLr::loo_roc_rate_calc] and defined as follows:\itemize{ 
#'          \item \eqn{Recall = TP / (TP + FN)} Recall is also know at the true positive rate.
#'          \item \eqn{Precision = TP / (TP + FP)}
#'          \item {Where, *TP* = number of true positive top assignments, 
#'          *FN* = number of false negative top assignments, 
#'          *FP* = number of false positive top assignments from the leave-one-out analysis}
#'          }
#' 
#' @returns a tibble with the minimum and maximum threshold for each repunit (aka reporting group)
#' 
#' @examples
#' sa_input <- rubias::self_assign(reference = GCLr::ex_baseline, gen_start_col = 5)
#' 
#' IA_thresholds(input = input, rates = NULL, min_rec = 0.80, min_pre = 0.95, min_thres = 0.50, thres_levels = seq(0.01, 0.99, by = 0.01))
#' 
#' @export
IA_thresholds <- function(sa_input, rates = NULL, min_rec = 0.80, min_pre = 0.95, min_thres = 0.50, thres_levels = seq(0.01, 0.99, by = 0.01)){
  
  group_names <- sa_input$repunit %>% 
    unique()
  
  if(is.null(rates)){
    
    rates <- GCLr::loo_rate_calc(data = sa_input, thres_levels = seq(0.01, 0.99, by = 0.01), group_names = group_names, ncores = parallel::detectCores())
    
  }
    
  output <- rates %>% 
      dplyr::filter(level >= min_thres, tpr >= min_rec, pre >= min_pre) %>% 
      dplyr::group_by(repunit) %>% 
      dplyr::summarize(min_threshold = min(level), max_threshold = max(level))
  
  return(output)

}