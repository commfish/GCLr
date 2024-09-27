#' Individual Principal Component Analysis (PCA)
#'
#' This function reads in a GENEPOP file and performs an individual Principal Component Analysis. 
#'
#' @param file the file path of a GENEPOP file, with .gen extension, containing individual genotypes. Make sure that each individual in the file has a unique identifier (e.g., popA_1, popA_2, popA_3). If using `GCLr::gcl2genepop()` to create the GENEPOP file, set `VialNums = TRUE`.
#' 
#' @param groupvec A numeric vector indicating the group affiliation of each population in `dist.mat`
#' 
#' @param groups A character vector of group names with a length equal to the maximum value in `groupvec`
#' 
#' @param ncomp The number of principal coordinates to include in the output (default = 3).
#' 
#' @details This function requires the adegenet and ade4 packages. 
#' For this function to work properly, the GENEPOP file must have a ".gen" extension (required by [adegenet::read.genepop()]), and 
#' each individual must have a unique identifier (see "file" argument description). The output of this function can be supplied to [GCLr::plot_indiv_pca()] to plot the first 3 principal components of the analysis in an interactive 3-dimensional scatter plot.
#' The variance explained by each principal coordinate are calculated from the eigen values in the raw PCA output (i.e., (pca.raw$eig)/sum(pca.raw$eig)). 
#' 
#' @seealso [ade4::dudi.pca()]
#' @seealso [GCLr::plot_indiv_pca()]
#' 
#' @return A named list with three elements:
#'  \itemize{
#'       \item \code{pca.raw}: the raw output from [ade4::dudi.pca()]
#'       \item \code{pca.df}: a tibble containing group, collection, indiv (unique individual identifier), and coordinates for each principal component.
#'       \item \code{var_exp}: the variance explained for all principal components (calculated from eigen values in pca.raw)
#'       }
#'
#' @examples
#'
#' group_info <- GCLr::ex_baseline %>%
#'   dplyr::mutate(collection = factor(collection, levels = unique(collection))) %>% 
#'   dplyr::group_by(repunit, collection) %>%
#'   dplyr::summarise(group = unique(repunit) %>% factor()) %>% 
#'   dplyr::arrange(collection)
#' 
#' file <- system.file("genepop", "ex_genepop.gen", package = "GCLr")
#' 
#' file <- "C:/Users/awbarclay/Documents/R/GCLr/inst/genepop/ex_genepop.gen" #DELETE for final version
#' 
#' groups <- group_info$group %>% unique()
#'
#' groupvec <- group_info$group %>%
#'   factor(levels = groups) %>%
#'   as.numeric()
#' 
#' indiv_pca(file = file, groupvec = groupvec, groups = groups, ncomp = 3)
#' 
#' @export
indiv_pca <- function(file, groupvec, groups, ncomp = 3){
  
  gp <- adegenet::read.genepop(file = file)
  
  scaled_freqs <- adegenet::scaleGen(gp, NA.method = "mean")
  
  ind_pca <- ade4::dudi.pca(scaled_freqs, center = FALSE, scale = FALSE, scannf = FALSE, nf = ncomp)
  
  group_info <- (gp@pop %>% janitor::tabyl()) %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(pop = gsub(pattern = "*_\\d+", replacement = "", x = .)) %>% 
    dplyr::mutate(group = groups[groupvec]) %>% 
    dplyr::select(pop, group)
  
  new_names <- tibble::as_tibble(ind_pca$li, rownames = "indiv") %>% 
    names() %>%  
    gsub(pattern = "Axis", replacement = "PC", x = .)
  
  ind_pca_df <- tibble::as_tibble(ind_pca$li, rownames = "indiv") %>% 
    purrr::set_names(new_names) %>% 
    dplyr::mutate(pop = gsub(pattern = "*_\\d+", replacement = "", x = indiv)) %>% 
    dplyr::left_join(group_info, by = "pop") %>% 
    dplyr::mutate(group = factor(group, levels = unique(group))) %>% 
    dplyr::select(group, collection = pop, indiv, dplyr::everything())#DF of coordinates
  
  var_exp <- (ind_pca$eig)/sum(ind_pca$eig) #Proportion of variance explained by each PC. Only the first three PCs will be used in plot
  
  output <- list(pca.raw = ind_pca, pca.df = ind_pca_df, var_exp = var_exp)
  
  return(output)  
  
}