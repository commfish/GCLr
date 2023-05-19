load_sillys <- function(path, sillyvec = NULL, rds = FALSE) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function loads .gcl objects saved with `dput`, it is a wrapper for `dget`.
  # It assumes all .txt files in `path` are .gcls and will load all of them,
  # unless you specify `sillyvec`, in which case it will just load those.
  # Also assumes files do not have ".gcl" in filename.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   path - character vector of where the .gcl objects you wish to load reside
  #   sillyvec - (optional) character vector of sillys in `path`
  #   rds - logical; if set to TRUE the function will load files in rds format.
  #                  if set to FALSE the function will load files produced by 'dput'
  #                  default is FALSE for backwards compatablity.
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a character vector of sillys loaded
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # setwd("V:/Analysis/1_SEAK/Sockeye/Mixture/Lynn Canal Inseason/2018/")
  # load_sillys(path = "Baseline genotpyes", sillyvec = c("SCKAT07E", "SCKAT07L"))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  if (!require("pacman")) {
    install.packages("pacman")
  }
  library(pacman)
  pacman::p_load(tidyverse)  #  Install packages, if not in library and then load them
  
  if (rds) {
    extension <- ".rds"
  } else {
    extension <- ".txt"
  }
  
  if (is.null(sillyvec)) {
    files_to_load <-
      base::list.files(
        path = path,
        pattern = extension,
        full.names = FALSE,
        recursive = FALSE
      )
    
  } else {
    files_to_load <- base::paste0(sillyvec, extension)
    
    if (!all(
      files_to_load %in% base::list.files(
        path = path,
        pattern = extension,
        full.names = FALSE,
        recursive = FALSE
      )
    )) {
      stop("Not all sillys in `sillyvec` are in `path`")
    }
  }
  
  if (length(files_to_load) == 0) {
    stop(base::paste0("There are no '", extension, "'files in the path provided"))
  }
  
  objects <- invisible(sapply(files_to_load, function(file) {
    obj <-
      stringr::str_replace(string = file,
                           pattern = extension,
                           replacement = ".gcl")
    
    if (rds) {
      base::assign(x = obj,
                   value = readRDS(file = paste(path, file, sep = "/")),
                   pos = 1)
    } else {
      base::assign(x = obj,
                   value = dget(file = paste(path, file, sep = "/")),
                   pos = 1)
    }
    
    obj
    
  }, USE.NAMES = FALSE))
  
  print(objects)
  
}