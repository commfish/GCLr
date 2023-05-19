add_tree_color <- function(tree, currentnames, treenames, groupvec, regioncol, regionpch = NULL, write_nexus = FALSE, file = NULL){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function takes a tree object produced by ape::nj() and adds group colors to the branch lengths for each reporting groups. 
  # Additionally, this function can used to modified the tip labels with new name and optional symbols.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   tree - a phylogenetic tree object (class = "phylo") produced by the ape package
  #
  #   currentnames - a character vector of the tip labels (pop names) in the tree object
  #   treenames - a character vector of new tip labels (pop names) to replace current names
  #   groupvec - numeric vector indicating the group affiliation for each population
  #   *Note: currentnames, treenames, and groupvec should be the same length and order.*
  #
  #   regioncol - is a vector of colors the same length as max(groupvec) (i.e. the number of groups).
  #               The colors can be hexadecimal, R color names, or the number corresponding to an R color name in colors() 
  #
  #   regionpch - is a vector with the same structure as "regioncol", but each element is an R pch (plot character) number. 
  #
  #   write_nexus - when TRUE the function calls on save_nexus_color_tree() to write out a Nexus tree file that can be opened 
  #                 in FigTree. FigTree can be downloaded from GitHub (https://github.com/rambaut/figtree/releases) 
  #                 *Note: Nexus files do not support R plot characters (pch), so they will not be included in the file.*
  #     
  #   file - when write_nexus = TRUE, the full file path for writing out a Nexus tree file including the .nex extension.
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a new tree object with colors and tip labels.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # attach("V:/Analysis/2_Central/Coho/Cook Inlet/2019/2019_Cook_Inlet_coho_baseline/2019_Cook_Inlet_coho_baseline.Rdata")
  # tree <- FstTree
  # currentnames <- tree$tip.label
  # treenames <- gsub(pattern = "\\)", x = gsub(pattern ="\\(", x = treenames104, replacement = ""), replacement = "")
  # groupvec <- groupvec_11
  # detach()
  # 
  # r_colors <- c("purple", "blue", "forestgreen","magenta", "red", "orange", "green", "yellow", "cyan", "gray", "brown") #R color names example
  # 
  # hex_colors <- c("#A020F0FF", "#0000FFFF", "#228B22FF", "#FF00FFFF", "#FF0000FF", "#FFA500FF", "#00FF00FF", "#FFFF00FF", "#00FFFFFF", "#BEBEBEFF", "#A52A2AFF") #Hexadecimal colors example
  # 
  # r_color_numbers <- match(r_colors, colors()) #R color numbers example
  # 
  # colortree <- add_tree_color(tree = tree,
  #               currentnames = currentnames,
  #               treenames = treenames,
  #               groupvec = groupvec,
  #               regioncol = r_color_numbers,
  #               regionpch = unique(groupvec),
  #               write_nexus = TRUE,
  #               file = "ColoredNexusTree.nex")
  # 
  # ape::plot.phylo(x = colortree$tree, edge.color = colortree$color, edge.width = 3, use.edge.length = T, show.tip.label = T, adj = .02, cex = .7, font = 1, label.offset = 0.002)
  # ape::tiplabels(pch = colortree$pch , offset = .0015, cex = .5)
  # 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(write_nexus == TRUE&is.null(file)){
    
    stop("A file path must be supplied when write_nexus = TRUE")
    
  }  
    
  if(write_nexus == TRUE&!is.null(regionpch)){
    
    warning("Nexus files do not support R plot characters (pch), so they were not included in the file.")
    
  }
    
  if(is.numeric(regioncol)){
    
    regioncol <- colors()[regioncol]
    
  }
  
  regioncol <- as.character(regioncol)
  
  if(grepl(pattern = "black", x = regioncol) %>% sum() > 0){
    
    warning("If you didn't know this already, using black as a regioncol will make the population branches look the same as the rest of the tree segments.")
    
  }
  
  xx <- tree
  pops1 <- currentnames
  pops2 <- treenames
  pops3 <- groupvec
  C <- length(pops1)

  xx$tip.label <- pops2[match(xx$tip.label, pops1)]

  regInd <- pops3[match(xx$tip.label, pops2)]
  cols <- regioncol[regInd]
  edgeind <- match(xx$edge[ , 2], 1:C)
  EdgeReg <- regInd[edgeind]


  all <- cbind(1:length(EdgeReg), xx$edge[ , 1], EdgeReg)

  edgereg <- EdgeReg
  MM <- length(EdgeReg)
  kk <- 1
  while(kk){
    
  NN <- dim(all)[1]
  
  if(NN==0) break()
  
    ind <- (apply(all[-NN, 2:3]==all[-1,2:3], 1, prod)==1)
    ind[is.na(ind)] <- F
    ind[1] <- F
    ind2 <- match(all[,1][ind],all[,1])
    all[ind2-1, 3] <- all[ind2, 3]
    edgereg[all[, 1][ind]-1] <- all[ind2, 3]
    all <- all[-c(ind2, ind2 + 1),] 
    
  }

  EdgeReg <- edgereg

  edgecol <- regioncol[EdgeReg]

  edgecol[is.na(edgecol)] <- "black"

  if(!is.null(regionpch)){
    
    mypch <- regionpch[regInd]
    
  } else{mypch <- NULL}
  
  if(write_nexus == TRUE ){
    
    save_nexus_color_tree(tree = tree, file = file, groupvec = groupvec, group.col = regioncol, tip.labs = xx$tip.label , sep = "_")
    
  }

  return(list(tree = xx, color = edgecol, pch = mypch))

}