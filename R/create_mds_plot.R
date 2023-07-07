#' Create an Interactive Multidimensional Scaling (MDS) Plot
#'
#' This function takes a pairwise matrix of genetic distances (e.g. Fst, CSE, Nei's) and creates an interactive MDS plot that can be save as an .html.
#'
#' @param file the file path (including html extension) to save plot as an html. If NULL, no file will be saved (defult = NULL)
#' 
#' @param dist.mat pairwise distance matrix
#' 
#' @param pop_names a character vector of population names corresponding to each population in `dist.mat` to add hover labels to the points on the plot, if no names are supplied the labels will default to dimnames(dist.mat)[[1]]
#' 
#' @param groupvec A numeric vector indicating the group affiliation of each population in `dist.mat`
#' 
#' @param group_names A character vector of group names with a length equal to the maximum value in `groupvec`
#' 
#' @param group_colors a character vector of colors corresponding to each group in `group_names`
#' 
#' @param title the title to place at the top of the plot; if NULL (default) no title will be added to the plot
#' 
#' @param labels weather to add labels the points on the MDS plot (default = FALSE)
#' 
#' @param showlegend whether to include a legend for the group colors (default = TRUE)
#' 
#' @return an interactive MDS plot 
#'
#' @examples
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' loci <- GCLr::ex_baseline[,-c(1:5)] %>%
#'   names() %>%
#'   gsub(pattern = "*\\.1", x = ., replacement = "") %>%
#'   unique()
#' 
#' dist.mat <- GCLr::pairwise_fst(sillyvec = sillyvec, loci = loci, inputfile = system.file("genepop", "ex_genepop.txt", package = "GCLr"), ncores = 4)
#' 
#' group_info <- GCLr::ex_baseline %>% 
#'   dplyr::group_by(collection) %>% 
#'   dplyr::summarise(group = unique(repunit) %>% factor())
#' 
#' GCLr::create_mds_plot(file = NULL, dist.mat = dist.mat, pop_names = group_info$collection, groupvec = group_info$group %>% as.numeric(), group_names = levels(group_info$group), group_colors = c("red", "blue", "green"), title = "Example MDS plot using pairwise Fst matrix", labels = TRUE, showlegend = TRUE)
#' 
#' @export
create_mds_plot <- function(file, dist.mat, pop_names, groupvec, group_names, group_colors, title = NULL, labels = FALSE, showlegend = TRUE){
  
  if(!length(unique(lengths(list(pop_names, groupvec, dimnames(dist.mat)[[1]])))) == 1){
    
    stop("pop_names, groupvec, and both dimensions of dist.mat must be the same length")
    
  }
  
  if(!length(unique(lengths(list(group_colors, group_names)))) == 1){
    
    stop("group_colors and group_names must be the same length")
    
  }
  
  if(!is.matrix(dist.mat)){
    
    stop("dist.mat must be a matrix object")
    
  }
  
  if(is.null(pop_names)){
    
    pop_names <- dimnames(dist.mat)[[1]]
    
  }
  
  xx <- cmdscale(dist.mat, k = 3)
  
  coords <- tibble::tibble(group_name = factor(group_names[groupvec], levels = group_names), PC1 = as.vector(xx[ , 1]), PC2 = as.vector(xx[ , 2]), PC3 = as.vector(xx[ , 3]), pop_names)
  
  p <- coords %>% 
    plotly::plot_ly(x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers", color = ~group_names[groupvec], colors = group_colors, text = ~pop_names, hoverinfo = "text") %>%   
    plotly::layout(title = title, showlegend = showlegend)
  
  if(!is.null(file)){htmlwidgets::saveWidget(plotly::as.widget(p), file)}
  
  return(p)
  
}