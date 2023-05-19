plot_quick_tree <- function(sillyvec, loci, inputfile, popnames = NULL, ncores = 4, groupvec = NULL, regioncol = NULL, cex = 1){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function creates a quick neighbor-joining Fst tree using the genepop and ape packages
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - a character vector of locus names
  # 
  #   inputfile - the file path of the genepop input file including .txt extension. If the file does not exist, then one will be created.
  #
  #   popnames - optional vector of population names corresponding to sillys in sillyvec to add as the dimnames of the output maxtrix. e.g., dimnames(ouput.matrix) <- list(popnames, popnames)
  #              If popnames is not supplied, sillyvec will be used as the dimnames e.g., dimnames(ouput.matrix) <- list(sillyvec, sillyvec)
  #   
  #   ncores - a numeric vector of length one indicating the number of cores to use
  #
  #   groupvec - numeric vector indicating the group affiliation for each population in sillyvec
  #
  #   regioncol - is a vector of colors the same length as max(groupvec) (i.e. the number of groups).
  #               The colors can be hexadecimal, R color names, or the number corresponding to an R color name in colors() 
  #
  #   cex - a numeric value giving the factor scaling of the tip labels (Character EXpansion).
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  #  A neighbor-joining phylogram based on pairwise Fst values.
  #
  #  Note: If groupvec is supplied and regioncol = NULL, groups will be colored using rainbow(n = max(groupvec))
  #        If groupvec and regioncol a both NULL, the tree will not be colored.
  #  
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   source(paste0(path.expand("~/R/"), "Functions.R"))#GCL functions
  #
  #   load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #
  #   plot_quick_tree(sillyvec = sillyvec31, loci = loci82, inputfile = "genepop/Susitna31pops82loci.txt", popnames = treenames31, ncores = 8, groupvec = groupvec4, regioncol = grcol4, cex = .75)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  pwfst <- pairwise_fst(sillyvec = sillyvec, loci = loci, inputfile = inputfile, popnames = popnames, ncores = ncores)
  
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