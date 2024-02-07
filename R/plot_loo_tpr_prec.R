#' Plot True Positive Rate Vs Precision from Leave-One-Out Analysis
#' 
#' This function takes the unmodified output from `GCLr::loo_roc_rate_calc()`, 
#' calculates precision, and plots the true positive rate (y-axis) by precision (x-axis)
#' in an interactive plot that can be used to determine assignment thresholds for 
#' individual assignment analyses.
#' 
#' @param input the unmodified output object from `GCLr::loo_roc_rate_calc()`
#' 
#' @param file the file name with .html extension for saving the plots.
#' 
#' @details -	TP = number of true positives, FN = number of false negatives, FP = number of false positives from the leave-one-out analysis.
#'            \itemize{
#'                  \item \eqn{TPR = TP / (TP + FN)}
#'                  \item \eqn{Precision = TP / (TP + FP)}
#'            }
#'            
#'          - Details about the leave-one-out analysis can be found on the \pkg{rubias} GitHub page: \url{https://github.com/eriqande/rubias} or
#'            in Moran and Anderson (2018).
#'             
#'          - citation: Moran, B.M. and E.C. Anderson. 2019. Bayesian inference from the conditional genetic stock 
#'                      identification model. Canadian Journal of Fisheries and Aquatic Sciences, 76(4): 551–560. 
#'                      
#' @returns This function returns an interactive plot of True positive rate (y-axis) and Precision (x-axis) faceted by reporting group. 
#'          Dashed red lines indicate the cutoffs for precision (>= 0.95) and true positive rate (>= 0.80). 
#'          The red point on the line indicates where the assignment threshold = 50%. If a file name is supplied, the function will save the plot to an .html file.
#'                      
#' @examples
#' group_names <- GCLr::ex_baseline$repunit %>% unique()
#' 
#' self_assign <- rubias::self_assign(reference = GCLr::ex_baseline, gen_start_col = 5)
#' 
#' input <- GCLr::loo_roc_rate_calc(data = self_assign, thres_levels = seq(0.01, .99, by = 0.005), group_names = group_names) # Note: make sure to include 0.50 in your threshold levels for 
#' 
#' plot_loo_tpr_prec(input, file = "~/example_plot.html") #Saves example plot to your documents folder
#'      
#' @export
plot_loo_tpr_prec <- function(input, file = NULL){
  
  repunits <- input$repunit %>% 
    unique()
  
  #look for the closest threshold to 50%
  fiftypct <- input %>% 
    dplyr::filter(level >= .5) %>% 
    dplyr::filter(level == min(level)) %>% 
    dplyr::pull(threshold) %>% 
    unique()
  
  dat <- input %>% 
    dplyr::mutate(tpr_marker = dplyr::case_when(threshold == fiftypct ~ tpr),
                  pre_marker = dplyr::case_when(threshold == fiftypct ~ pre),
                  threshold_marker = dplyr::case_when(threshold == fiftypct ~ threshold))
  
  # Get minimum precision and true positive rate for setting limits
  min_prec <- dat$pre %>% min() %>% round(digits = 3)
  
  min_tpr <- dat$tpr %>% min() %>% round(digits = 3)
  
  point_dat <- dat %>% 
    dplyr::mutate(repunit = factor(repunit, levels = repunits)) %>% 
    dplyr::filter(threshold_marker == threshold)
  
  ggp <- dat %>% 
    dplyr::mutate(repunit = factor(repunit, levels = repunits)) %>% 
    ggplot2::ggplot(ggplot2::aes(x = pre, y = tpr, label = threshold)) + 
    ggplot2::geom_path(color = "blue") + 
    ggplot2::geom_point(data = point_dat , ggplot2::aes( x = pre, y = tpr), color = "red")+
    ggplot2::xlab("Precision")+
    ggplot2::ylab("True Positive Rate")+
    ggplot2::geom_vline(xintercept = 0.95, color = "red", linetype = 2)+
    ggplot2::geom_hline(yintercept = 0.80, color = "red", linetype = 2)+
    ggplot2::ylim(c(min_tpr, 1.02))+
    ggplot2::xlim(c(min_prec, 1.02))+
    ggplot2::facet_wrap(~repunit, nrow = (length(repunits)/3) %>% ceiling(), ncol = 3)
  

  p <- plotly::ggplotly(ggp)
  

  
 
 if(!is.null(file)){
   
   htmlwidgets::saveWidget(p, file = file)
   
 }
 
  return(p)
  
}