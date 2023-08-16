#' Get metadata from GT-seq projects
#' 
#' This function pulls a GT-seq project metadata report from Loki
#' 
#' @param project_name a character vector of one or more GT-seq project names as they are spelled in LOKI. e.g., c("P014","P015","P016") (see details)
#' 
#' @param dir the directory for writing output file with "\\" separator between folders. e.g., "V:\\Analysis\\Staff\\Andy Barclay\\R"
#' 
#' @param file_name the name of the csv file. (see details)
#' 
#' @param username your user name for accessing LOKI through R, a.k.a your state username
#' 
#' @param password our password for accessing LOKI through R.
#' 
#' @param open.file logical,if set to TRUE the CSV file will open after it has been written. Default is FALSE.
#'
#' @returns a tibble of metadata by locus (crosstab/wide format) and writes out to one or more .csv files (see details)
#' 
#' @details
#' This function is only intended for pulling metadata for GT-seq projects. If a project name is supplied that is not a GT-seq project, the function will still produce files, but there will be NA's in most of the fields.
#' Reports that exceed the maximum number of rows for Excel (1048575 rows max) will be split into multiple files with a number added to the file name (e.g. test_metadat_1.csv, test_metadat_2.csv)
#'
#' @examples
#' \dontrun{
#'   metadata <- get_gtseq_metadata(project_name = "K158", dir = "C:\\Users\\awbarclay\\Documents", file_name = "test_metadat", username = "awbarclay", password = "************", open.file = FALSE)
#' }
#' 
#' @export
get_gtseq_metadata <- function(project_name, dir, file_name, username, password, open.file = FALSE){  
 
  # Recording function start time
  start.time <- Sys.time()
  
  # Remove extension if people add it to file_name
  if(length(grep("*.csv", file_name, value = FALSE)) == 1){
    
    file_name <- stringr::str_split(string = file_name, pattern = ".csv")[[1]][1]
    
  }
  
  # Setting default java.parameters
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  # Creating java query
  gnoqry <- paste0("SELECT LAB_PROJECT_ID, LAB_PROJECT_NAME, SILLY_CODE, FK_FISH_ID, LOCUS, POSITIONS, HAPLO_ALLELES, HAPLO_COUNTS, GENOTYPE, SNP_ALLELES, PROBES, FWD_PRIMER_SEQ, DNA_PLATE_ID, DNA_PLATE_WELL_POS FROM AKFINADM.R_GTSEQ_GENO_METADATA GENO WHERE LAB_PROJECT_NAME IN (", paste0("'", project_name, "'", collapse = ","), ")") #The view for this query (AKFINADM.R_GTSEQ_GENO_METADATA) was created for this function on 5/12/20.
 
  # Pull genotypes and concatenate alleles into one column with "/" separator
  dataAll0 <- RJDBC::dbGetQuery(con, gnoqry) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(LAB_PROJECT_NAME == project_name)

  # Disconnect from LOKI
  discon <- RJDBC::dbDisconnect(con)
  
  # Put in crosstab format
  dataAll <- dataAll0 %>% 
    dplyr::arrange(LOCUS, SILLY_CODE, FK_FISH_ID)

  # Checking to see if no data was pulled for a given project. If true, the function will stop and print an error message in the console.
  if(nrow(dataAll)==0 & exists("project_name")){
    
    stop("No data exist for the supplied lab project name. Check to make sure the name is spelled correctly. Also, project names for genotypes imported prior to October 2016 may not be included in LOKI lookup table.")
    
  }
  
  # Write data to an excel csv file with no more than 1,048,575 rows
  loc <- dataAll %>% 
    dplyr::pull(LOCUS) %>% 
    unique()
  
  n_indiv <- dim(dataAll %>% dplyr::filter(LOCUS == loc[1]))[1]
  
  xls <- 1048575 #for making sure final file is no longer than the Excel row limit
  
  nloci <- floor(xls/n_indiv)#Maximum number of loci to include in each file
  
  rows <- nloci*n_indiv #Max number of rows per file so loci aren't split among multiple files
    
  end <- length(dataAll$FK_FISH_ID)
  
  N <- base::ceiling(end/rows)#Number of files to write
  
  for (i in 1:N){
    
    start <- rows*i-rows+1
    
    stop <- if(i*rows > end){end} else(i*rows)
    
    readr::write_csv(dataAll[start:stop, ], file = paste(dir,"/", file_name, "_", i, ".csv", sep=''))
    
  }

  # Open CSV file in excel if open.file is set to TRUE
  if(open.file == TRUE){
    
    shell(paste(dir, "/", file_name, "_", i, ".csv", sep=''), wait = FALSE)
    
    }
  
  # Calculate the time it took the function to pull and write data and print time to console
  stop.time <- Sys.time() 
  
  message(paste0("CSV file writen to:", tools::file_path_as_absolute(dir)))
  
  print(stop.time-start.time)
  
  # Return data to assign to an object
  return(dataAll)
  
}

