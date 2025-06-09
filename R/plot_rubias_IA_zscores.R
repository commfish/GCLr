#' Plot distribution of `rubias` individual assignment z-scores.
#' 
#' This function reads in a `rubias` individual posteriors csv file produced by [GCLr::run_rubias_mix()] and produces a histogram of the individual assignment (IA) z-scores.
#' 
#' @param mixnames A character vector of mixture names to include in the IA summary.
#' @param path Character vector of where to find output from each mixture as a .csv (created by [GCLr::run_rubias_mix()]).
#' @param bins number of histogram bins (default = 30)
#' @param out_file Character vector of the file path to save the histograms as a pdf file including pdf extension. If NULL (default), no pdf will be produced.
#' 
#' @return a histogram of z-scores for each mixture in `mixnames`
#'   
#' @examples
#' \dontrun{
#' path <- "V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Mixture/UCI_sockeye_2020_postseason/rubias/output"
#' mixnames <- c("DriftDW_20", "DriftCorr_20", "UpperSub_20", "Kasilof600ft_20", "WestKalgin_20", "Eastern_20", "General_North_20", "General_South_20")
#' out_file <- "~/my_indiv_assign_zscore_plots.pdf"
#' 
#' GCLr::plot_rubias_IA_zscores(path = path, mixnames = mixnames, out_file = out_file)
#' }
#' @export
plot_rubias_IA_zscores <- function(path, mixnames, bins = 30, out_file = NULL){
  
  mixnames <- ind_assign$mixture_collection %>% 
    unique()
  
  plots <- pbapply::pblapply(mixnames, function(mix){
    
    mix.assign <- readr::read_csv(file = paste0(path, "/", mix, "_indiv_posteriors.csv"), show_col_types = FALSE)
    
    mix.assign %>% 
      dplyr::group_by(indiv) %>% 
      dplyr::slice_max(z_score, n = 1) %>% 
      ggplot2::ggplot(ggplot2::aes(x = z_score)) + 
      ggplot2::geom_histogram(bins = bins) + 
      ggplot2::theme_bw()+
      ggplot2::ggtitle(label = "Distribution of indiviudal assignment z-scores",
                       subtitle = paste0("mixture = ", mix)) 
    
  }) %>% purrr::set_names(mixnames)
  
  if(!is.null(out_file)){
    
    print(paste0("Saving a PDF of z-score histograms to: '", normalizePath(out_file), "'"))
    
    pdf(file = out_file)
    
    plots
    
    dev.off()
    
  }
  
  plots
  
}