#' Split a Combined Chip Import File into Multiple Files
#' 
#' This function takes a combined Biomark import file and splits it up into  
#' multiple files so genotypes can be imported into Loki in smaller batches (see details)
#' 
#' @param files the file path to the Biomark files including .csv extension
#' @param nlines the maximum number of lines you want in each file (default = 9216)
#' 
#' @details
#' The function reads in a Biomark file and splits it up into multiple files with a maximum of number of lines = `nlines`
#' and writes them out in the same folder as the original file with letters at the end of the file names. For example,
#' if you have a file named "S254 03_04 Combined Chip Run.csv" with 18,432 (192 samps X 24 loci X 4 chips) lines (not including the header lines) and you want to break it up into files with 9,216 lines each (96 samps X 96 loci), the function will produce
#' 2 files, and name them "S254 03_04 Combined Chip Run_a.csv" and "S254 03_04 Combined Chip Run_b.csv"
#' 
#' 
#' @examples
#' #' \dontrun{
#' files <- list.files(path = "C:/Users/awbarclay/Documents/biomark split tests", pattern = "*Combined Chip Run*", full.names = TRUE)
#' 
#' split_combchip_loki_import(files = files, nlines = 9216)
#' 
#' }
#' 
#' @export
split_combchip_loki_import <- function(files, nlines = 9216){
  
  lapply(files, function(file){
    
    my.file <- suppressWarnings(readr::read_csv(file, col_names = FALSE, skip_empty_rows = FALSE, col_types = readr::cols(.default = "c")))
    
    header <- my.file[1:15, ]
    
    fields <- my.file[16, ]
    
    dat <- my.file[-c(1:16), ]
    
    last.line <- dim(dat)[1]
    
    nfiles <- ceiling(last.line/nlines)
    
    start <- 1
    
    end <- nlines
    
    file_name <- strsplit(file, split = ".csv")[[1]]
    
    for(i in 1:nfiles){
      
      outfile <- dplyr::bind_rows(header, fields, dat[start:end, ])
      
      readr::write_csv(outfile, paste0(file_name, "_", letters[i], ".csv"), col_names = FALSE, na = "", eol = "\r\n") 
      
      start <- start+nlines
      
      if(i == nfiles-1){
        
        end <- dim(my.file)[1]
        
      }else{end <- end + nlines}
      
    }
    
  })
  
}