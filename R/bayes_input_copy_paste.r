#' Copy and Paste BAYES Input Files
#'
#' This function copies the BAYES input files (.ctl, .mix) from the source directory
#' and pastes them into the target directory to run BAYES.
#'
#' @param origindir Original directory of files.
#'   A character vector indicating where BAYES files are located.
#' @param targetdir Target directory structure.
#'   A character vector indicating where to paste BAYES files (i.e., on the server).
#' @param mixvec Mixtures to be moved.
#'   A character vector indicating which mixtures to move.
#' @param i Bayes folder identifier.
#'   Indicates which "Bayes" folder to start putting files into.
#'   Default is to put the first mixture into "Bayes A", second into "Bayes B", etc.
#'   This function assumes that the Bayes folders on the V drive are listed alphabetically.
#'
#' @return This function has no explicit return value. It copies and pastes the files onto the server.
#'
#' @Author: Kyle Shedd
#' 
#' @examples
#' \dontrun{
#' bayes_input_copy_paste(origindir = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/BAYES/Late August 89loci",
#'                        targetdir = "C:/Users/krshedd/BAYES", 
#'                        mixvec = dget(file = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects/LateAugustMixtures2014Strata.txt"),
#'                        i = "A")
#'
#' }
#'
#' @export
#' 
bayes_input_copy_paste <- function(origindir, targetdir, mixvec, i = "A"){
  
  controlfiles <- list.files(path = paste(origindir, "/Control", sep = ""), pattern = ".ctl", full.names = TRUE, recursive = FALSE) # All control files in origindir
 
  mixturefiles <- list.files(path = paste(origindir, "/Mixture", sep = ""), pattern = ".mix", full.names = TRUE, recursive = FALSE) # All mixture files in origindir

  filestocopy = sapply(mixvec, function(mix) {
    
    c(controlfiles[grep(pattern = mix, x = controlfiles)], mixturefiles[grep(pattern = mix, x = mixturefiles)])
    
    }, simplify=FALSE) # Just the control/mixture files for mixvec
  
  j <- which(LETTERS == i) - 1
  
  invisible(sapply(seq(mixvec), function(mix) {
    
    file.copy(from = filestocopy[[mix]], to = paste(targetdir, "/Bayes ", LETTERS[mix + j], sep = ""))
    
    }))
  
}