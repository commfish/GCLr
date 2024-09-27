#' Plot a Principal Component Analysis (PCA)
#'
#' This function takes the raw output from [GCLr::indiv_pca] and plots the coordinates for the first 3 principal components in an interactive 3-dimensional scatter plot. 
#'
#' @param pca the raw output from [GCLr::indiv_pca]
#' 
#' @param group_colors a character vector of colors corresponding to each group in the raw PCA output (groups = pca$pca.df$group %>% unique())
#' 
#' @param title the title to place at the top of the plot; if NULL (default) no title will be added to the plot
#' 
#' @return A 3-dimensional scatter plot of PCs 1-3
#' 
#' @seealso [GCLr::indiv_pca()]
#'
#' @examples
#' \dontrun{
#' # Paste this code into your console to run this example.
#' group_info <- GCLr::ex_baseline %>%
#'   dplyr::mutate(collection = factor(collection, levels = unique(collection))) %>% 
#'   dplyr::group_by(repunit, collection) %>%
#'   dplyr::summarise(group = unique(repunit) %>% factor()) %>% 
#'   dplyr::arrange(collection)
#' 
#' file <- system.file("genepop", "ex_genepop.gen", package = "GCLr")
#' 
#' groups <- group_info$group %>% unique()
#'
#' groupvec <- group_info$group %>%
#'   factor(levels = groups) %>%
#'   as.numeric()
#' 
#' pca_out <- indiv_pca(file = file, groupvec = groupvec, groups = groups, ncomp = 3)
#' 
#' plot_indiv_pca(pca = pca_out, group_color = c("blue", "green", "red"), title = NULL)
#' 
#' }
#' 
#' @export
plot_indiv_pca <- function(pca, group_color, title = NULL){
  
  #Check to make sure the pca input contains at least 3 PCs
  npc <- pca$pca.df %>% 
    names() %>% 
    grep(pattern = "PC", x = .) %>% 
    length()
  
  if(npc < 3){
    
    stop("The PCA input must contain coordinates for at least 3 PCs.")
    
  }
  
  groups = pca$pca.df$group %>% 
    unique()
  
  
  group_color <- (group_color %>% 
                    GCLr::col2hex())$hex %>% 
    purrr::set_names(groups)
  
  plot_df <- pca$pca.df %>% 
    dplyr::mutate(group_col = group_color[group])
 
  p <- plot_df %>% 
    plotly::plot_ly(x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers", color = ~group, colors = group_color, text = ~indiv, hoverinfo = "text") %>% 
    plotly::layout(title = title, 
                   showlegend = TRUE, 
                   scene = list(
                     xaxis = list(title = paste0("PC1(", round(pca$var_exp[1]*100, 1), "%)")),
                     yaxis = list(title = paste0("PC2(", round(pca$var_exp[2]*100, 1), "%)")),
                     zaxis = list(title = paste0("PC3(", round(pca$var_exp[3]*100, 1), "%)"))))
  
  return(p)  
  
}