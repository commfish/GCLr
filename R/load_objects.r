#' Load R Objects
#'
#' This function loads R objects saved with [dput()] or [saveRDS()], serving as a wrapper for [dget()] and [readRDS()].
#'
#' @param path A character vector specifying the path where the objects to load reside.
#' @param pattern Optional character vector argument to manually specify a pattern (i.e., specific object or regular expression). Accepts a list of patterns if you want to load multiple objects (e.g., \code{c("pattern1", "pattern2")}).
#' @param rds Logical; if set to TRUE, the function will load files in RDS format. If set to FALSE, the function will load text files produced by [dput()]. The default is FALSE for backward compatibility.
#' 
#' @return A character vector of the objects loaded into the workspace.
#' 
#' @examples
#' \dontrun{
#' 
#' GCLr::load_objects(path = "V:/Analysis/1_SEAK/Sockeye/Mixture/Lynn Canal Inseason/2018/Objects", pattern = "^loci") # Loads loci objects from "Objects" directory
#' 
#' }
#' 
#' @export
load_objects <- function(path, pattern = NULL, rds = FALSE) {
  
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