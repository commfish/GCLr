save_nexus_color_tree <- function(tree, file, groupvec, group.col, tip.labs, sep = "_"){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function is a wrapper for ape::write.nexus() and creates a Nexus format tree file including tip labels and branch colors, 
  #   which can be opened and modified in FigTree.
  #
  #   FigTree can be downloaded from GitHub: https://github.com/rambaut/figtree/releases
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   tree - a phylogenetic tree object (class = "phylo") produced by the ape package
  #
  #   file - the file path with file name including .nex extension
  #
  #   groupvec - a numeric vector indicating the group affiliation of each tip label
  #
  #   group.col - character vector of R colors() length of max(groupvec) to color the branches
  #
  #   tip.labs - a character vector same length as groupvec to label the tips of the tree.
  #              Note: tip.labs cannot have spaces to work in FigTree.  If tip.labs contains spaces,
  #                    this function will replace them with an underscore.
  #
  #   sep - delimiter used to replace spaces in the tip.labs; if spaces are left in the tip.labs (called TAXALABLES in the NEXUS file), when FigTree opens the file 
  #         an error message will pop up saying the number of taxa is greater than ntaxa.
  #              
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #  A nexus tree file
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #
  #   save_nexus_color_tree(tree = FstTree, file = "Susitna_Chinook_tree.nex", groupvec = groupvec, group.col = grcol, tip.labs = Final_Pops$location)
  # 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  if(!length(tree$tip.label)==length(tip.labs)){
    
    stop(paste0("The tip.labs supplied is not the same length as the number of tips on the tree: ", length(tree$tip.label)))
    
  }
  
  if(!length(groupvec)==length(tree$tip.label)){
    
    stop(paste0("The groupvec supplied is not the same length as the number of tips on the tree: ", length(tree$tip.label)))
    
  }
  
  if(!length(group.col)==max(groupvec)){
    
    stop(paste0("The colors supplied are not the same length as the number of groups in groupvec: ", max(groupvec)))
    
  }
  
  tip.labs <- gsub(x = tip.labs, pattern = " ", replacement = sep)
  
  tree$tip.label <- tip.labs

  ape::write.nexus(phy = tree, file = file)

  #Function to get Hexadecimal color codes
  GetColorHexadecimal <- function(color){
  
    c <- col2rgb(color)
    
    cbind(color, sprintf("#%02X%02X%02X", c[1], c[2], c[3]), sprintf("%3d %3d %3d", c[1], c[2], c[3]))
    
  }
    
  if(is.numeric(group.col)){
  
    group.col <- colors()[group.col]

  }
          
  GrCol <- group.col
          
  black <- GetColorHexadecimal("black")[2]
          
  hexcol <- sapply(GrCol, function(col){
    
    GetColorHexadecimal(col)[2]
    
    })
       
  rle <- rle(groupvec)

  PopCol <- unlist(sapply(seq(length(rle$values)), function(col){
    
    rep(hexcol[rle$values[col]], rle$lengths[col])
    
    })) #Vector of colors the same order and length as groupvec
       
  names(PopCol) <- seq(length(tip.labs))
       
  nexfile <- scan(file, what = '', sep = "\n", blank.lines.skip = FALSE)
       
  treerow <- grep(pattern = "END;", nexfile)[2]-1
       
  tree <- nexfile[treerow]
       
  for(pop in 1:length(tip.labs)){

    tree <- sub(x = tree, 
                pattern = paste("[(]", pop, ":", sep = ''),
                replacement = paste("(", tip.labs[pop], "[&!color=", PopCol[pop], "]:", sep = ''))
      
    tree <- sub(x = tree, pattern = paste("[,]", pop, ":", sep = ''),
                replacement = paste(",", tip.labs[pop], "[&!color=", PopCol[pop], "]:", sep = ''))
      
  }
      
  tree <- gsub(x=tree,pattern="):",replacement=paste(")[&!color=",black,"]:",sep=''))
      
  nexfile[treerow] <- tree
      
  tstart <- grep(x = nexfile, pattern = "	TRANSLATE")
      
  tend <- tstart+length(tip.labs)+1

  nexfile <- nexfile[-c(tstart:tend)]
      
  write.table(nexfile, file = file, quote = FALSE, row.names = FALSE, col.names = FALSE )
 
}

