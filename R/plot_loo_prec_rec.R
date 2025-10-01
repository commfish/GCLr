#' Produce a precision-recall curve plot
#' 
#' This function takes the unmodified output from [rubias::self_assign], plots the scaled likelihood for each baseline sample (y-axis) by reporting group or population in a box plot (optional), and plots recall (aka *true positive rate*) (x-axis) by precision (y-axis)
#' in an interactive plot that can be used to assess reporting groups and assignment thresholds for individual assignment analyses.
#' 
#' @param sa_input the unmodified output object from [rubias::self_assign]
#' 
#' @param rates optional; the unmodified output object from [GCLr::loo_rate_calc] created using the supplied sa_input. If left NULL, the rates will be calculated from sa_input. Supplying a rates object will make this function run a lot faster.
#' 
#' @param group_colors a vector of colors the same length as the number of groups in the input file (i.e., `input$repunit %>% unique()`).
#' 
#' @param by_group logical; whether you want the box plot by reporting group (TRUE) or population (FALSE); default = TRUE
#' 
#' @param min_rec the recall (true positive rate) cutoff for adding a vertical dashed line to the precision plot; default = 0.80
#' 
#' @param min_pre the precision cutoff for adding a horizontal dashed line to the precision plot; default = 0.95
#' 
#' @param min_thres minimum assignment likelihood threshold for adding a horizontal line at on the boxplot; default = 0.50
#' 
#' @param thres_levels a vector of assignment thresholds to include on precision plot (only used when `rates = NULL`). 
#' 
#' @param plots the plots do you want the function to produce. There are three choices:
#'              \itemize{
#'                \item{\emph{facet}: }{An iteractive precision-recall curve plot with all reporting groups (i.e., faceted plot}
#'                \item{\emph{paired}: }{Interative plots for each reporting group containing a box plot of assignement likelihoods on the left hand side and a precision-recall curve plot on the right hand side.}
#'                \item{\emph{both}: }{Facet plot first, then paired plots.}
#'              }
#' 
#' @param ncores the number of cores to use. Default is to use all available cores.
#' 
#' @details Details about the leave-one-out analysis (loo) can be found on the \pkg{rubias} GitHub page: \url{https://github.com/eriqande/rubias} or
#'          in Moran and Anderson (2018). Citation: Moran, B.M. and E.C. Anderson. 2019. Bayesian inference from the conditional genetic stock 
#'          identification model. Canadian Journal of Fisheries and Aquatic Sciences, 76(4): 551–560.
#'          Recall and precision are are calculated from the loo results with [GCLr::loo_rate_calc] and defined as follows:\itemize{ 
#'          \item \eqn{Recall = TP / (TP + FN)} Recall is also know at the true positive rate.
#'          \item \eqn{Precision = TP / (TP + FP)}
#'          \item {Where, *TP* = number of true positive top assignments, 
#'          *FN* = number of false negative top assignments, 
#'          *FP* = number of false positive top assignments from the leave-one-out analysis}
#'          }
#'            
#' @returns When `plot = "facet"`, the function returns an interactive line plot of True positive rate (y-axis) and Precision (x-axis) faceted by reporting group (aka repunit), with dashed red lines indicating the cutoffs for precision (>= 0.95) and true positive rate (>= 0.80). 
#'          The red point on the line indicates where the assignment threshold = 50%.
#'          When `plot = "paired"`, the function returns a list of interactive plots for each reporting group in the leave-on-out analysis. 
#'          The left hand side is a box plot of the scaled likelihood of each baseline sample assigning to the reporting group, with a dashed red line indicating the 50% assignment threshold. 
#'          The right hand side is a line plot of True positive rate (y-axis) and Precision (x-axis), with dashed red lines indicating the cutoffs for precision (>= 0.95) and true positive rate (>= 0.80). 
#'          The red point on the line indicates where the assignment threshold = 50%.
#'          When `plot = "both"`, the function returns a list of plots with the faceted TPR vs Precision plot as the first element and the paired plots for each reporting group for remaining elements.
#'          When no rates are supplied, an object named 'rates' is assigned to your global environment.
#'              
#' @note The plots produced by this function are not intended for publications; however, the code within this function can be copied and modified to produce plots tailored for publications.            
#'   
#' @seealso [GCLr::loo_rate_calc]     
#'                
#' @examples
#' sa_input <- rubias::self_assign(reference = GCLr::ex_baseline, gen_start_col = 5)
#' 
#' # Plots are returned as a list when more than one plot is produced 
#' plots <- GCLr::plot_loo_prec_rec(sa_input, group_colors = c("red", "green", "blue"), by_group = FALSE, min_rec = 0.80, min_pre = 0.95, min_thres = 0.50, thres_levels = seq(0.01, 0.99, by = 0.01), plots = c("facet", "paired", "both")[3], ncores = parallel::detectCores())
#' 
#' plots[[1]] # This is how you would look at the first plot in the list
#'      
#' @export
plot_loo_prec_rec <- function(sa_input, rates = NULL, group_colors, by_group = TRUE, min_rec = 0.80, min_pre = 0.95, min_thres = 0.50, thres_levels = seq(0.01, 0.99, by = 0.01), plots = c("facet", "paired", "both")[3], ncores = parallel::detectCores()){
  
  repunits <- sa_input$repunit %>% 
    unique()
  
  if(length(group_colors) != length(repunits)){
    
    stop("The group_colors supplied to not match the number of groups in sa_input.")
    
  }
  
  group_colors <- purrr::set_names(group_colors, repunits) # Add group names to color vector
  
  #Get pop info and arrange by repunit so all pops from a repunit are side by side. 
  #Pop numbers will be preserved by the order in which they appeared in sa_input.

  pop_info <- sa_input %>% 
    dplyr::select(collection, repunit) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(pop_order = seq_along(collection), repunit = factor(repunit, levels = repunits)) %>% 
    dplyr::arrange(repunit)
  
  groupvec <- pop_info$repunit
  
  popvec <- pop_info$collection #Vector of pop names
  
  pop_col <- group_colors[groupvec] %>% 
    purrr::set_names(popvec)#color for each pop
  
  groupvec_rle <- rle(as.character(groupvec))$lengths %>% #This object is used to create the segmented line of group colors below the box plot
    purrr::set_names(rle(as.character(groupvec))$values)
  
  #Box plot data
  
  box_dat <- sa_input %>% 
      dplyr::group_by(indiv, collection, repunit, inferred_repunit) %>%
      dplyr::summarise(repu_scaled_like = sum(scaled_likelihood), .groups = "drop_last")
  
  #Precision plot data
  
  if(is.null(rates)){ #Calculating rates if none supplied.
    
    rates <- GCLr::loo_roc_rate_calc(data = sa_input, thres_levels = thres_levels, group_names = repunits, ncores = ncores)
     
    assign(x = "rates", value <- rates, envir = .GlobalEnv) # Assign the rates object to global env.
    message("'rates' object assinged to global environment")
    
  }
  
  pre_dat0 <- dplyr::mutate(rates, repunit = factor(repunit, levels = repunits))

  fiftypct <- pre_dat0 %>% #filter for the closest threshold to min_thres
    dplyr::mutate(level = as.character(level)) %>% 
    dplyr::mutate(level = as.numeric(level)) %>% # For some reason filter wasn't working for certain min_thres values. Converting to character and then back to numeric solves the issue. WTF?!
    dplyr::filter(level >= min_thres) %>% 
    dplyr::filter(level == min(level)) %>% 
    dplyr::pull(threshold) %>% 
    unique()
  
  pre_dat <- pre_dat0 %>% 
    dplyr::mutate(tpr_marker = dplyr::case_when(threshold == fiftypct ~ tpr),
                  pre_marker = dplyr::case_when(threshold == fiftypct ~ pre),
                  threshold_marker = dplyr::case_when(threshold == fiftypct ~ threshold)) %>% 
    dplyr::rename(precision = pre, recall = tpr)
  
  point_dat <- pre_dat %>% #For point on line at the 50% threshold
    dplyr::mutate(repunit = factor(repunit, levels = repunits)) %>% 
    dplyr::filter(threshold_marker == threshold)
  
  # Faceted plots
  if(plots %in% c("facet", "both")){
    
    plot1 <- plotly::ggplotly(pre_dat %>% 
      dplyr::mutate(repunit = factor(repunit, levels = repunits)) %>% 
      ggplot2::ggplot(ggplot2::aes(x = recall, y = precision, label = threshold, color = repunit)) + 
      ggplot2::geom_path() + 
      ggplot2::scale_color_manual(values = group_colors) +
      ggplot2::geom_point(data = point_dat , ggplot2::aes( x = recall, y = precision), color = "red")+
      ggplot2::ylab("Precision")+
      ggplot2::xlab("Recall")+
      ggplot2::geom_hline(yintercept = min_pre, color = "red", linetype = 2)+
      ggplot2::geom_vline(xintercept = min_rec, color = "red", linetype = 2)+
      ggplot2::facet_wrap(~repunit)) %>% 
      plotly::layout(showlegend = FALSE)
    
  } 
  
  if(plots %in% c("paired", "both")){
    
    # Create plots for each repunit
    plot2 <- lapply(repunits, function(repu){ 
      
      # Box plot
      plot_data <- box_dat %>% 
        dplyr::mutate(Population = factor(collection, levels = popvec)) %>% 
        dplyr::filter(inferred_repunit == repu) %>%
        dplyr::mutate(repunit = factor(repunit, levels = repunits), inferred_repunit = factor(inferred_repunit, levels = repunits)) %>% 
        dplyr::arrange(collection) %>% 
        dplyr::mutate(`Scaled Likelihood` = repu_scaled_like, group_col = group_colors[repunit]) 
      
      hline <- function(y = 0, color = "red"){ #This function adds a horizontal line to the plot
        
        list(type = "line", x0 = 0, x1 = .9, xref = "paper", y0 = y, y1 = y, line = list(color = color, dash = "dash"))
        
      }
      
      if(by_group == FALSE){ #By pop
        
        p1 <- plotly::plot_ly(plot_data, x = ~Population, y = ~repu_scaled_like, colors = group_colors) %>% 
          plotly::add_boxplot(color = ~repunit, marker = list(line = list(color = "black", width = .5))) %>% 
          plotly::layout(showlegend = FALSE, shapes = hline(0.5), yaxis = list(title = "Scaled Likelihood"), xaxis = list(showticklabels = FALSE) )
        
      } else{ #By group
        
        p1 <- plotly::plot_ly(plot_data, x = ~repunit, y = ~repu_scaled_like, colors = group_colors) %>%
          plotly::add_boxplot(color = ~repunit, marker = list(line = list(color = "black", width = .5))) %>% 
          plotly::layout(showlegend = FALSE, shapes = hline(0.5), yaxis = list(title = "Scaled Likelihood"), xaxis = list(title = "Reporting Group") )
        
      }
      
      #Precision-recall plot
      p2 <- plotly::ggplotly(pre_dat %>% 
                               dplyr::filter(repunit == repu) %>% 
                               ggplot2::ggplot(ggplot2::aes(x = recall, y = precision, label = threshold)) + 
                               ggplot2::geom_path(color = group_colors[repu]) + 
                               ggplot2::geom_point(data = point_dat %>% dplyr::filter(repunit == repu), ggplot2::aes( x = recall, y = precision), color = "red")+
                               ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 0, vjust = .1), 
                                     axis.text.y = ggplot2::element_text(size = 12),
                                     axis.title = ggplot2::element_text(size = 14),
                                     title = ggplot2::element_text(size = 12, face = "bold"),
                                     legend.position = "none")+
                               ggplot2::ylab("Precision")+
                               ggplot2::xlab("Recall")+
                               ggplot2::geom_hline(yintercept = min_pre, color = "red", linetype = 2)+
                               ggplot2::geom_vline(xintercept = min_rec, color = "red", linetype = 2)+
                               ggplot2::ylim(c(0,1)) +
                               ggplot2::ggtitle(paste0("Assignment to: ", repu)))
        
      
     suppressWarnings(plotly::subplot(p1, p2, shareY = FALSE, shareX = FALSE, titleX = TRUE, titleY = TRUE, widths = c(0.80, 0.20)) %>% 
       plotly::layout(yaxis = list(tickfont = list(size = 10)), 
                      yaxis2 = list(tickfont = list(size = 10)), 
                      xaxis = list(domain = c(0, 0.70), tickfont = list(size = 10)), 
                      xaxis2 = list(domain = c(0.8, 1), tickfont = list(size = 10))))# Combine box and precision plots
      
    })
    
  }

  if(plots == "facet"){
      return(plot1)
  }
  
  if(plots == "paired"){
      return(plot2)
  }
 
  if(plots == "both"){
     return(list(plot1, plot2))
  }
  
}