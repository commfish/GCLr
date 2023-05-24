#' Plot Quick Phylogenetic Tree
#'
#' This function quickly creates a neighbor-joining Fst tree using the \pkg{genepop} and \pkg{ape} packages.
#'
#' @param sillyvec a vector of silly codes without the ".gcl" extension.
#' 
#' @param loci a vector of locus names.
#' 
#' @param inputfile the file path of the genepop input file including the .txt extension. If the file does not exist, then one will be created.
#' 
#' @param popnames an optional vector of population names corresponding to sillyvec in order to add them as the dimnames of the output matrix. If not supplied, sillyvec will be used as the dimnames.
#' 
#' @param ncores a numeric value indicating the number of cores to use
#' 
#' @param groupvec a numeric vector indicating the group affiliation for each population in sillyvec
#' 
#' @param regioncol a vector of colors of the same length as max(groupvec) (i.e., the number of groups). The colors can be hexadecimal, R color names, or the number corresponding to an R color name in \code{colors()}.
#' 
#' @param cex a numeric value giving the factor scaling of the tip labels (Character Expansion).
#'
#' @return a neighbor-joining phylogram based on pairwise Fst values.
#'
#' @examples
#' source(paste0(path.expand("~/R/"), "Functions.R"))  # GCL functions
#' load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
#' plot_quick_tree(sillyvec = sillyvec31, loci = loci82, inputfile = "genepop/Susitna31pops82loci.txt", popnames = treenames31, ncores = 8, groupvec = groupvec4, regioncol = grcol4, cex = 0.75)
#'
#' @import ape
#' 
#' @export

plot_quick_tree <- function(sillyvec, loci, inputfile, popnames = NULL, ncores = 4, groupvec = NULL, regioncol = NULL, cex = 1){
  
  pwfst <- GCLr::pairwise_fst(sillyvec = sillyvec, loci = loci, inputfile = inputfile, popnames = popnames, ncores = ncores)
  
  tree <- ape::nj(pwfst)
  
  if(is.null(groupvec) & is.null(regioncol)){
    
    regioncol <- "black"
    
  }
  
  if(!is.null(groupvec) & is.null(regioncol)){
    
    regioncol <- rainbow(n = max(groupvec))
    
  }
  
  if(is.null(groupvec) & !is.null(regioncol)){
    
    stop("If regioncol is supplied, then groupvec must be supplied to color the tree.")
    
  }
  
  colortree <- add_tree_color(tree = tree, currentnames = dimnames(pwfst)[[1]], treenames = popnames, groupvec = groupvec, regioncol = regioncol)
  
  plot_color_tree(color.tree = colortree, cex = cex)
  
}