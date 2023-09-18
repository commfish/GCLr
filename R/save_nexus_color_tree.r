#' Save Nexus File that includes a Phylogenetic Tree. 
#'
#' This function is a wrapper for [ape::write.nexus()] and creates a Nexus format tree file including tip labels and branch colors, which can be opened and modified in FigTree. FigTree can be downloaded from GitHub: [FigTree Releases](https://github.com/rambaut/figtree/releases)
#'
#' @param tree a phylogenetic tree object (`class = "phylo"`) produced by the \pkg{ape} package
#' @param file the file path with file name including .nex extension
#' @param groupvec a numeric vector indicating the group affiliation of each tip label
#' @param group.col a character vector of R [grDevices::colors()] with a length of `max(groupvec)` to color the branches
#' @param tip.labs A character vector with the same length as groupvec to label the tips of the tree. tip.labs cannot have spaces to work in FigTree. If tip.labs contains spaces, this function will replace them with an underscore.
#' @param sep Delimiter used to replace spaces in tip.labs. If spaces are left in tip.labs (called TAXALABLES in the NEXUS file) when FigTree opens the file, an error message will pop up saying the number of taxa is greater than ntaxa.
#'
#' @examples
#' 
#' groupvec <- GCLr::ex_baseline %>%
#'   dplyr::group_by(collection) %>%
#'   dplyr::filter(dplyr::row_number()==1) %>%
#'   dplyr::pull(repunit) %>%
#'   factor() %>%
#'   as.numeric()
#' 
#' sillyvec <- GCLr::base2gcl(ex_baseline)
#' 
#' loci <- GCLr::ex_baseline[,-c(1:5)] %>%
#'   names() %>%
#'   gsub(pattern = "*\\.1", x = ., replacement = "") %>%
#'   unique()
#' 
#' pwfst <- GCLr::pairwise_fst(sillyvec = sillyvec, loci = loci, inputfile = system.file("genepop", "ex_genepop.txt", package = "GCLr"), ncores = parallel::detectCores())
#' 
#' FstTree <- ape::nj(pwfst)
#' 
#' GCLr::save_nexus_color_tree(tree = FstTree, file = path.expand("~/example_tree.nex"), groupvec = groupvec, group.col = c("green", "red", "blue"), tip.labs = paste0("pop", 1:14))
#' 
#' @export
save_nexus_color_tree <- function(tree, file, groupvec, group.col, tip.labs, sep = "_"){

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

  if(is.numeric(group.col)){
  
    group.col <- colors()[group.col]

  }
          
  GrCol <- group.col
          
  black <- GCLr::col2hex("black")$hex
          
  hexcol <- sapply(GrCol, function(col){
    
    GCLr::col2hex(col)$hex
    
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
      
  tree <- gsub(x = tree, pattern = "):", replacement = paste(")[&!color=",black,"]:", sep = ''))
      
  nexfile[treerow] <- tree
      
  tstart <- grep(x = nexfile, pattern = "	TRANSLATE")
      
  tend <- tstart+length(tip.labs)+1

  nexfile <- nexfile[-c(tstart:tend)]
      
  write.table(nexfile, file = file, quote = FALSE, row.names = FALSE, col.names = FALSE )
 
}