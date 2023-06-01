#' Plot a Colored Phylogenetic Tree
#'
#' This function is a wrapper for [ape::plot.phylo()]. It takes a tree object created by [GCLr::add_tree_color()] and plots it with an optional scale.
#'
#' @param color.tree a color tree object created with [GCLr::add_tree_color()].
#' 
#' @param rm.neg.branch logical; if TRUE (default), negative branch lengths are removed.
#' 
#' @param type a character string specifying the type of phylogeny to be drawn. It must be one of "phylogram" (the default), "cladogram", "fan", "unrooted", "radial", or any unambiguous abbreviation of these.
#' 
#' @param scale logical; whether to add a scale to the tree plot.
#' 
#' @param show.tip.label logical; if TRUE (default), labels are placed at the tip of each branch.
#' 
#' @param edge.width a numeric vector giving the width of the branches of the plotted phylogeny. If fewer widths are given than the length of edge, then these are recycled.
#' 
#' @param adj a numeric specifying the justification of the text strings of the labels: 0 (left-justification), 0.5 (centering), or 1 (right-justification). This option has no effect if type = "unrooted". If NULL (the default), the value is set with respect to the direction (see details in [ape::plot.phylo()] documentation).
#'
#' @param cex a numeric value giving the factor scaling of the tip labels (Character Expansion).
#' 
#' @param font an integer specifying the type of font for the labels: 1 (plain text), 2 (bold), 3 (italic, the default), or 4 (bold italic)
#' 
#' @param label.offset a numeric giving the space between the nodes and the tips of the phylogeny and their corresponding labels. This option has no effect if type = "unrooted".
#' 
#' @param gen.dist the genetic distance used to produce the tree. If "fst", "FST", or "Fst" is supplied, the "F" will be italicized and the "ST" will be subscripted.
#' 
#' @param ... further arguments to be passed to [ape::plot.phylo()]
#'
#' @return This function plots a phylogenetic tree.
#'
#' @examples
#' load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
#' Fst36 <- read_genepop_pwfst(file = "GENEPOP/Sustina36pops84loci.MIG", popnames = Pooling_1_pops$location)
#' FstTree <- ape::nj(Fst36)
#' groupvec <- factor(Pooling_1_pops$group, unique(Pooling_1_pops$group)) %>% as.numeric()
#' colortree <- add_tree_color(tree = FstTree, currentnames = FstTree$tip.label, treenames = FstTree$tip.label, groupvec = groupvec, regioncol = grcol)
#' plot_color_tree(color.tree = colortree, rm.neg.branch = TRUE, type = "phylogram", scale = TRUE, show.tip.label = TRUE, edge.width = 3, adj = 0.02, cex = 0.7, font = 1, label.offset = 0.001, gen.dist = "FST")
#'
#' @import ape
#' 
#' @aliases PlotColorTree.GCL
#' 
#' @export

plot_color_tree <- function(color.tree, rm.neg.branch = TRUE, type = "phylogram", scale = TRUE, show.tip.label = TRUE, edge.width = 3, adj = NULL, cex = 1, font = 1, label.offset = 0, gen.dist = "FST", ...){

  if(rm.neg.branch == TRUE){
    
    color.tree$tree$edge.length = pmax(0, color.tree$tree$edge.length) #Get rid of negative branches  
    
  }
  
  ape::plot.phylo(x = color.tree$tree, type = type, edge.color = color.tree$color, edge.width = edge.width, use.edge.length = TRUE, show.tip.label = show.tip.label, adj = adj, cex = cex, font = font, label.offset = label.offset, ...)
  
  if(scale == TRUE){
    
    axis(1)  #Adds scale to bottom of plot
    
    if(gen.dist %in% c("FST", "fst", "Fst")){
      
      scale.lab <- expression(italic(F)[ST])
      
    }else{
      
      scale.lab <- gen.dist
      
    }
    
    mtext(text = scale.lab, side = 1, cex = 1.5, outer = F, line = 3)
    
  }
    
}
#' @rdname plot_color_tree
#' @export
PlotColorTree.GCL <- plot_color_tree  