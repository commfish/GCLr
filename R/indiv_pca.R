#' Individual Principal Component Analysis (PCA)
#'
#' This function reads in a GENEPOP file and performs an individual Principal Component Analysis. 
#'
#' @param file the file path of a GENEPOP file with .gen extension containing individual genotypes. Make sure that each individual in the file has a unique identifier (e.g., popA_1, popA_2, popA_3). If using [GCLr::gcl2genepop()] to create the GENEPOP file, set `VialNums = TRUE`.
#' 
#' @param groupvec A numeric vector indicating the group affiliation of each collection in the GENEPOP file.
#' 
#' @param groups A character vector of group names with a length equal to the maximum value in `groupvec`
#' 
#' @param scale_freqs A logical vector of length 1 that controls whether individual allele frequencies are scaled (i.e., variance = 1); default = `FALSE`. `TRUE` gives you DAPC (discriminant analysis of principal components), `FALSE` gives you covariance PCA.
#' 
#' @param ncomp The number of principal coordinates to include in the output (default = 3).
#' 
#' @details This function requires the adegenet and ade4 packages. 
#' For this function to work properly, the GENEPOP file must have a ".gen" extension (required by [adegenet::read.genepop()]), and 
#' each individual must have a unique identifier (see "file" argument description). The output of this function can be supplied to [GCLr::plot_indiv_pca()] to plot the first 3 principal components of the analysis in an interactive 3-dimensional scatter plot.
#' The variance explained by each principal coordinate are calculated from the eigen values in the raw PCA output (i.e., (pca.raw$eig)/sum(pca.raw$eig)).
#' When `scale_freqs` is `TRUE` you get DAPC (discriminant analysis of principal components, i.e., how can pre-defined groups best be separated?), which will upweight any rare alleles and can cause individuals with rare alleles to show up as outliers.
#' When `scale_freqs` is `FALSE`, you get normal, covariance PCA (i.e., what are the major axes of overall genetic variation), which will not produce extreme outliers if individuals have rare alleles.
#' It is recommended to use covariance PCA to explore population structure and DAPC when groups are known and you want maximal discrimination.
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
#' groups <- group_info$group %>% unique()
#'
#' groupvec <- group_info$group %>%
#'   factor(levels = groups) %>%
#'   as.numeric()
#' 
#' indiv_pca(file = file, groupvec = groupvec, groups = groups, scale_freqs = FALSE, ncomp = 3)
#' 
#' @export
indiv_pca <- function(file, groupvec, groups, scale_freqs = FALSE, ncomp = 3){
  
  gp <- adegenet::read.genepop(file = file)
  
  if(scale_freqs) {

    scaled_freqs <- adegenet::scaleGen(gp, NA.method = "mean")  # scales and centers allele frequencies (i.e., mean = 0, variance = 1)
    
    ind_pca <- ade4::dudi.pca(scaled_freqs, center = FALSE, scale = FALSE, scannf = FALSE, nf = ncomp)  # DAPC; alleles contribute equally, used for uSATs
    
  } else {
    
    raw_freqs <- adegenet::tab(gp, freq=TRUE, NA.method="mean")  # raw allele counts
    
    ind_pca <- ade4::dudi.pca(raw_freqs, center=TRUE, scale=FALSE, scannf = FALSE, nf = ncomp)  # covariance PCA based on allele frequencies, rare alleles downweighted; centers each allele (i.e., mean = 0)
    
  }
  
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