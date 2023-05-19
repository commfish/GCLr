bayes_input_copy_paste <- function(origindir, targetdir, mixvec, i = "A"){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function will copy the BAYES input files (.ctl, .mix) from the V drive
  # and paste them onto the server to run BAYES.
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # origindir = original directory of files
  #   ~ A character vector indicating where BAYES files are on the V drive.
  # targetdir = target directory structure
  #   ~ A character vector indicating where you want BAYES files to go (i.e. on server).
  # mixvec = mixtures you want to move
  #   ~ A character vector indicating which mixtures you want to move.
  # i = "A"
  #   ~ Indicates which "Bayes" folder to start with putting files into.
  #     Default is to put 1st mixture into "Bayes A", second into "Bayes B", etc.
  #     This function assumes you have your Bayes folders on the V drive listed
  #     alphabetically.
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # No output, just copies files and pastes them on to the server for you.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   source("https://raw.githubusercontent.com/krshedd/GCL-R-Scripts/master/bayes_input_copy_paste.R")
  #   
  #   bayes_input_copy_paste(origindir="V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/BAYES/Late August 89loci",
  #                           targetdir="C:/Users/krshedd/BAYES", 
  #                           mixvec=dget(file="V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects/LateAugustMixtures2014Strata.txt"),
  #                           i = "A")
  #    
  # Created by Kyle Shedd Fri Mar 27 16:05:57 2015
  # Updated by Kyle Shedd Tue Apr 05 15:37:24 2016 in order to generalize and better document
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
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