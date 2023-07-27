#' Load Silly Files into the Global Environment
#'
#' This function loads silly files into the global environment. Silly files can
#' either be in '.txt' or '.rds' format. The function lists all the files in the
#' given \code{path} with the specified \code{extension}, loads them into the 
#' global environment as '.gcl' objects, and returns the names of the loaded objects.
#'
#' @param path A character string specifying the directory path where the silly 
#' files are located.
#'
#' @param sillyvec A character vector of specific silly file names (without the 
#' '.txt' or '.rds' extension) to be loaded. If provided, only the files listed 
#' in \code{sillyvec} will be loaded from the \code{path} directory.
#'
#' @param rds Logical. If TRUE, the function expects '.rds' files, otherwise, it 
#' assumes '.txt' files. Default is FALSE.
#'
#' @return A character vector containing the names of the loaded '.gcl' objects 
#' in the global environment.
#'
#'
#' @details The function uses the 'pacman' package to ensure that required 
#' packages are installed and loaded. If 'pacman' is not installed, it will 
#' automatically install it.
#'
#' The function loads the silly files into the global environment as '.gcl' 
#' objects. If the \code{rds} argument is set to TRUE, 'readRDS' from the 
#' 'base' package will be used to read the '.rds' files. Otherwise, the 'dget' 
#' function from the 'base' package will be used to read the '.txt' files.
#'
#' If the \code{sillyvec} argument is provided, the function will only load the 
#' specified files from the \code{path} directory. If any of the files in 
#' \code{sillyvec} are not present in the \code{path} directory, an error will be
#' thrown.
#'
#' If no files with the specified \code{extension} are found in the \code{path} 
#' directory, the function will raise an error indicating that there are no 
#' '.txt' or '.rds' files in the provided path.
#'
#' @examples
#' \dontrun{
#' # Load all '.txt' files from the 'silly_files' directory
#' load_sillys("silly_files")
#'
#' # Load specific '.rds' files from the 'data' directory
#' load_sillys("data", sillyvec = c("file1", "file2"), rds = TRUE)
#' }
#'
#' @export

load_sillys <- function(path, sillyvec = NULL, rds = FALSE) {  
  
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