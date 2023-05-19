load_objects <- function(path, pattern = NULL, rds = FALSE) {
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function loads R objects saved with `dput` or 'saveRDS', this is 
  # a wrapper for 'dget' and 'readRDS'.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   path - character vector of where the objects you wish to load reside
  #
  #   pattern - optional argument so you can manually specify a pattern (i.e., specific object)
  #             accepts a list of patterns if you want multiple objects: c("pattern1", "pattern2")
  #
  #   rds - logical; if set to TRUE the function will load files in rds format.
  #                  if set to FALSE the function will load test files produced by 'dput'
  #                  default is FALSE for backwards compatablity.
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a character vector of objects loaded
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # setwd("V:/Analysis/1_SEAK/Sockeye/Mixture/Lynn Canal Inseason/2018/")
  # load_objects(path = "Objects", pattern = "^loci") - just loads loci objects from "Objects" dir
  # load_objects(path = "Objects", pattern = c("^loci", "sample_size") - loads loci and sample size objects from "Objects" dir
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (!require("pacman")) {install.packages("pacman")}
  library(pacman)
  pacman::p_load(tidyverse)  #  Install packages, if not in library and then load them
  
  if (rds) {
    extension <- ".rds"
  } else {
    extension <- ".txt"
  }
  
  if (is.null(pattern)) {
    files_to_load <-
      list.files(
        path = path,
        pattern = paste0("*", extension),
        full.names = FALSE
      )
    
    if (length(files_to_load) == 0) {
      stop(paste0("There are no '", extension, "'files in the path provided"))
    }
    
  } else {
    files_to_load <-
      list.files(
        path = path,
        pattern = paste0("*", extension),
        full.names = FALSE
      )
    
    # files_to_load <-  stringr::str_subset(files_to_load, pattern = pattern)
    files_to_load <-
      stringr::str_subset(files_to_load, pattern = paste(pattern, collapse = "|")) #  multiple pattern matches
    
    if (length(files_to_load) == 0) {
      stop(
        paste0(
          "There are no '",
          extension,
          "'files containing the pattern ",
          "'",
          pattern,
          "' in the path provided"
        )
      )
    }
    
  }
  
  if (length(files_to_load) == 0) {
    stop(paste0("No files contain the pattern ", "'", pattern, "'"))
  }
  
  objects <- invisible(sapply(files_to_load, function(file) {
    obj <- unlist(strsplit(x = file, split = extension))
    
    if (rds) {
      assign(x = obj,
             value = readRDS(file = paste(path, file, sep = "/")),
             pos = 1)
      
    } else {
      assign(x = obj,
             value = dget(file = paste(path, file, sep = "/")),
             pos = 1)
      
    }
    
    obj
    
  }, USE.NAMES = FALSE))
  
  print(objects)
  
}