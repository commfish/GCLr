#' Plot Baseline Evaluation Results
#'
#' This function takes the baseline evaluation summary object produced by [GCLr::summarize_rubias_baseline_eval()] and produces plots of correct allocation for each test_group with summary statistics placed in the upper left corner of each plot.
#'
#' @param summary A baseline evaluation summary list produced by [GCLr::summarize_rubias_baseline_eval()].
#' @param file A tibble produced by [GCLr::base_eval_sample_sizes()] containing the following variables: 
#'    \itemize{
#'        \item \code{test_group} 
#'        \item \code{scenario}
#'        \item \code{repunit}
#'        \item  \code{samps}
#'    }
#' @param method A character string indicating the rubias output to summarize. Select one of three choices: 
#'    \itemize{
#'       \item {"MCMC"(plot MCMC output)}
#'       \item{"PB" (plot bias corrected output)}
#'       \item{"both" (plot both outputs)}
#'       }
#' @param test_groups A character vector of group names to include in the plots. This also sets the order in which they are plotted. If \code{test_groups} is not supplied, all test groups in \code{summary} will be plotted (see details).
#' @param group_colors A character vector of R \code{colors()} the same length as \code{test_groups}. If \code{group_colors} is not supplied, colors will be automatically selected using \code{rainbow()}.
#' 
#' @return If \code{method == "MCMC"} or \code{"PB"}, the function produces a single page pdf file containing the plots and returns the faceted plots. If \code{method == "both"}, the function produces a two-page pdf file containing the plots for both methods and no plots are returned.
#'
#' @examples
#' load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
#' require(tidyverse)
#' tests <- sample_sizes %>% group_by(test_group, scenario) %>% summarize(test_group = test_group %>% unique(), scenario = scenario %>% unique(), .groups = "drop_last")#Total of 510 tests  #
#' mixvec <- tests %>% unite(col = "mixvec", test_group, scenario, sep ="_" ) %>% pull()
#' path <-  "V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/rubias/3groups/output"
#' summary <- summarize_rubias_base_eval (mixvec = mixvec, sample_sizes = sample_sizes, method = "both", group_names = NULL, group_names_new = NULL, groupvec = NULL, groupvec_new = NULL, path = path, alpha = 0.1, burn_in = 5000, threshold = 5e-7, ncores = 8)
#' 
#' plot_baseline_eval(summary = summary, file = "Baseline_eval_plots.pdf", method = "PB", test_groups = groups3, group_colors = c("green", "magenta", "red"))
#'
#' @details The solid line on the plots indicates where the true proportion equals the estimated proportion. The two dotted lines indicate where the estimates fall within +/- 10% of the true proportion.
#' If a large number of \code{test_groups} are supplied, the faceted plots will become too small to see on one page. If so, you may need to supply a subset of \code{test_groups} to plot.
#' This function is not intended for producing publication-ready plots; however, the code from this function can be copied and modified to produce plots formatted for publication.
#' 
#' @import magrittr
#' @import ggplot2
#' @import dplyr
#' 
#' @aliases plot_baseline_eval_summary.GCL
#' 
#' @export
plot_baseline_eval <- function(summary, file, method = c("MCMC", "PB", "both"), test_groups = NULL, group_colors = NULL){
  
  # Check methods
  method_check <- summary$estimates$method %>% 
    unique() %>% 
    sort()
  
  if(method=="both"&!sum(method_check==sort(c("PB", "MCMC")))==2){
    
    stop("The baseline evalution summary supplied does not contain restults for both methods.")
    
  }
  
  if(method=="MCMC"&!sum(method_check%in%"MCMC")==1){
    
    stop("The baseline evalution summary supplied does not contain MCMC restults.")
    
  }
  
  if(method=="PB"&!sum(method_check%in%"PB")==1){
    
    stop("The baseline evalution summary supplied does not contain bias corrected (PB) results restults.")
    
  }
  
  # Get test_groups if none are supplied
  if(is.null(test_groups)){
    
    test_groups <- summary$summary_stats$test_group %>% 
      unique()
    
    }
  
  # Check test_groups
  test_group_check <- summary$estimates$test_group %>% 
    unique() %>% 
    as.character()
  
  if(!sum(test_group_check%in%test_groups)==length(test_groups)){
    
    stop("The baseline evaluation summary supplied does not contain results for all test_groups")
    
  }
 
  # Get rainbow of colors if no group_colors are supplied
  if(is.null(group_colors)){
    
    group_colors <- rainbow(n = length(test_groups))
    
  }
  
  # Color check
  if(!length(group_colors)==length(test_groups)){
    
    stop(paste0("The number of group_colors supplies is not the same length as test_groups: ", length(test_groups)))
    
  }
  
  pdf(file, height = 9, width=13)   
  
  # Both methods
  if(method=="both"){
    
    sapply(c("MCMC", "PB"), function(meth){
      
      plot <- summary$estimates %>%
        dplyr::filter(method==meth, test_group==repunit) %>% 
        dplyr::filter(test_group %in% test_groups) %>% 
        dplyr::left_join(summary$summary_stats, by = c("test_group", "method")) %>%
        dplyr::mutate(test_group = factor(test_group, levels = unique(as.character(test_group))),
                      repunit = factor(repunit, levels = unique(as.character(repunit)))) %>% 
        ggplot2::ggplot(ggplot2::aes(x = true_proportion, y = mean, colour = repunit)) +
        ggplot2::geom_point() +
        ggplot2::geom_linerange(ggplot2::aes(ymin = lo5CI, ymax = hi95CI))+
        ggplot2::geom_abline(intercept = 0, slope = 1) +
        ggplot2::geom_abline(intercept = 0.1, slope = 1, lty = 2) +
        ggplot2::geom_abline(intercept = -0.1, slope = 1, lty = 2) +
        ggplot2::scale_colour_manual(name = "Reporting Group", values = group_colors) +
        ggplot2::geom_text(ggplot2::aes(x = .3, y = 1, label = paste0("RMSE:", round(RMSE, digits = 3))), color = "black", size = 3)+ 
        ggplot2::geom_text(ggplot2::aes(x = .3, y = .94, label = paste0("Bias:", round(Mean_Bias, digits = 3))), color="black", size = 3)+
        ggplot2::geom_text(ggplot2::aes(x = .3, y = .88, label = paste0("90% Within: ", round(100*`90%_within`, 1), "%")), color = "black", size = 3)+
        ggplot2::geom_text(ggplot2::aes(x = .3, y = .82, label = paste0("Within Interval: ", round(100*Within_Interval, 1), "%")), color = "black", size = 3)+
        ggplot2::facet_wrap(~ test_group) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none", strip.text.x = ggplot2::element_text(size = 16), panel.spacing.y = ggplot2::unit(3, "lines"))+
        ggplot2::xlab("True Proportion") +
        ggplot2::ylab("Posterior Mean Reporting Group Proportion") +
        ggplot2::ggtitle(paste0("Baseline evaluation test results: ", meth), subtitle = paste0(summary$estimates$total_samps %>% unique() %>% ceiling(), " sample mixtures")) 
      
      print(plot)
      
    })
    
  } else{
    
    # One method
    meth <- method
    
    plot <- summary$estimates %>%
      dplyr::filter(method==meth, test_group==repunit) %>% 
      dplyr::filter(test_group %in% test_groups) %>% 
      dplyr::left_join(summary$summary_stats, by = c("test_group", "method")) %>%
      dplyr::mutate(test_group = factor(test_group, levels = unique(as.character(test_group))),
                    repunit = factor(repunit, levels = unique(as.character(repunit)))) %>% 
      ggplot2::ggplot(ggplot2::aes(x = true_proportion, y = mean, colour = test_group)) +
      ggplot2::geom_point() +
      ggplot2::geom_linerange(ggplot2::aes(ymin = lo5CI, ymax = hi95CI))+
      ggplot2::geom_abline(intercept = 0, slope = 1) +
      ggplot2::geom_abline(intercept = 0.1, slope = 1, lty = 2) +
      ggplot2::geom_abline(intercept = -0.1, slope = 1, lty = 2) +
      ggplot2::scale_colour_manual(name = "Reporting Group", values = group_colors) +
      ggplot2::geom_text(ggplot2::aes(x = .3, y = 1, label = paste0("RMSE:", round(RMSE, digits = 3))), color = "black", size = 3)+ 
      ggplot2::geom_text(ggplot2::aes(x = .3, y = .94, label = paste0("Bias:", round(Mean_Bias, digits = 3))), color="black", size = 3)+
      ggplot2::geom_text(ggplot2::aes(x = .3, y = .88, label = paste0("90% Within: ", round(100*`90%_within`, 1), "%")), color = "black", size = 3)+
      ggplot2::geom_text(ggplot2::aes(x = .3, y = .82, label = paste0("Within Interval: ", round(100*Within_Interval, 1), "%")), color = "black", size = 3)+
      ggplot2::facet_wrap(~ test_group) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none", strip.text.x = ggplot2::element_text(size = 16), panel.spacing.y = ggplot2::unit(3, "lines"))+
      ggplot2::xlab("True Proportion") +
      ggplot2::ylab("Posterior Mean Reporting Group Proportion") +
      ggplot2::ggtitle(paste0("Baseline evaluation test results: ", meth), subtitle = paste0(summary$estimates$total_samps %>% unique() %>% ceiling(), " sample mixtures"))
    
    print(plot)
    
  }
  
  dev.off()
  
  if(!method=="both"){plot}
  
}

#' @rdname plot_baseline_eval
#' @export
plot_baseline_eval_summary.GCL <- plot_baseline_eval