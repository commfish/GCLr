#' Copy and Paste BAYES Output Files
#'
#' This function copies the BAYES output files (.BO1, .BOT, .RGN, .SUM, .CLS) from the server
#' and pastes them onto the target directory (i.e., V drive). It also removes all BAYES input
#' files from the server. The function recursively searches the original directory for BAYES files
#' that match the mixvec.
#'
#' @param origindir Original directory of files.
#'   A character vector indicating where BAYES output files are located on the server.
#' @param targetdir Target directory structure.
#'   A character vector indicating where to paste BAYES files (i.e., V drive).
#' @param mixvec Mixtures to be moved.
#'   A character vector indicating which mixtures to move.
#'
#' @return This function has no explicit return value. It cuts and pastes BAYES output files
#'   from the server to the target directory, and deletes the BAYES input files from the server.
#'
#' @examples
#' \dontrun{
#' bayes_output_copy_paste(origindir = "C:/Users/krshedd/BAYES",
#'                         targetdir = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/BAYES/Late August 89loci",
#'                         mixvec = dget(file = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects/LateAugustMixtures2014Strata.txt"))
#'
#' }
#'
#' @export
#' 
bayes_output_copy_paste <- function(origindir, targetdir, mixvec){
  
  # Removing input files (".ctl" and ".mix")
  controlfiles = list.files(path = origindir, pattern = ".ctl", full.names = TRUE, recursive = TRUE) # All control files in origindir
  
  mixturefiles = list.files(path = origindir, pattern = ".mix", full.names = TRUE, recursive = TRUE) # All mixture files in origindir
  
  filestoremove = sapply(mixvec, function(mix) {
    
    c(controlfiles[grep(pattern = mix, x = controlfiles)], mixturefiles[grep(pattern = mix, x = mixturefiles)])
    
    }, simplify = FALSE) # Just the control/mixture files for mixvec 
  
  invisible(lapply(filestoremove, file.remove)) # Remove input files from server
  
  # Copy/paste output files (".BO1", ".BOT", ".RGN", ".SUM", and ".CLS")
  outputfiles = unlist(sapply(c(".BO1", ".BOT", ".RGN", ".SUM", ".CLS"), function(outfile) {list.files(path = origindir, pattern = outfile, full.names = TRUE, recursive = TRUE)})) # All output files in origindir
  
  filestocopy = sapply(mixvec, function(mix) {
    
    outputfiles[grep(pattern = mix, x = outputfiles)]
    
    }, simplify = FALSE) # Create list of output files by mixvec
  
  names(filestocopy) <- mixvec
  
  if (!all(mixvec %in% names(filestocopy))) {
    
    stop("'mixvec' object does not match all output files!!!\nThis 'stop' prevents you from accidently deleting these files, because they will not be copied properly!!!\nIf mixvec is a named vector, then the vector is used, not the names!!!")
    
    }
  
  invisible(sapply(mixvec, function(mix) {file.copy(from = filestocopy[[mix]], to = paste(targetdir, "/Output/", mix, sep = ""))})) # Move to appropriate V drive directory
  
  invisible(lapply(filestocopy, file.remove)) # Remove output files from server
  
}