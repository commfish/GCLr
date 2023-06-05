#' Get ASL Data From Loki
#'
#' This function connects to Loki and pulls the fish information (ASL data) for each silly in sillyvec. This pulls the same information as the "ASL Import" report in OceanAK.
#'
#' @param sillyvec a character vector of silly codes
#' 
#' @param username your state user name
#' 
#' @param password your password used to access Loki. Contact Eric Lardizabal if you don't have a password for Loki.
#' 
#' @param file The file path, including .csv extension, for writing out a csv file of the output. If set to NULL, no file will be saved.
#' 
#' @param import.vars Do you only want the variables allowed for the ASL importer? (default: TRUE)
#'
#' @return If `import.vars = TRUE`, the function outputs a tibble containing the following variables:
#'    \itemize{
#'        \item \code{Collection ID}
#'        \item \code{Fish ID}
#'        \item \code{Freshwater Age}
#'        \item \code{Ocean Age}
#'        \item \code{Sex}
#'        \item \code{Length}
#'        \item \code{Weight}
#'        \item \code{Scale Card Number}
#'        \item \code{Scale Card Position}
#'        \item \code{ASL Number}
#'    }
#' If `import.vars = FALSE`, a variable for `silly code` is added to the output tibble.
#'
#' @examples
#' get_asl_data(sillyvec = c("PPORTCH21", "PBARABCR21"), username = "awbarclay", password = password, import.vars = FALSE, file = "C:/Users/awbarclay/Documents/R/test_fish_table.csv")
#'
#' @details
#' The output of this function can be used for creating an import file for the Loki ASL Data Importer. The the importer requires that all columns are in the correct order and spelled correctly. The importer will not work if you include silly code in the import file.  When writing the import file, make sure it does not contain NAs When writing the import file, make sure it does not contain NAs. (see example below) # Also, when using readr::write_csv to write out the tissue import file, make sure to change the eol argument to from the default \verb{\\n} to \verb{\\r\\n} or the importer
#' Also, if you use [readr::write_csv()] to write out the tissue import file, make sure to change the eol argument to from the default \verb{\n} to \verb{\r\n} or the importer will give you an error message about the header names.  e.g.,  \verb{import_file %>% write_csv(file = "ImportFile.csv", na = "", eol = "\r\n")} 
#' This function requires an OJDBC driver object, which is an object in the GCLr package called [GCLr::drv]. 
#'
#' @aliases ASL_Import.GCL
#' 
#' @export 

get_asl_data <- function(sillyvec, username, password, file = NULL, import.vars = TRUE){
 
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  qry <- paste("SELECT * FROM AKFINADM.V_GEN_SAMPLED_FISH_TISSUE WHERE SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")", sep = "") #Query the tissues table to get collection IDs. Need to see if Eric can add silly code to the fish table.
  
  silly_collection_ID <- RJDBC::dbGetQuery(con, qry) %>% 
    dplyr::select(SILLY_CODE, FK_COLLECTION_ID) %>% 
    dplyr::distinct()#This gets the collection IDs.
  
  collection_IDs <- silly_collection_ID$FK_COLLECTION_ID
  
  qry <- paste0("SELECT * FROM AKFINADM.GEN_SAMPLED_FISH WHERE FK_COLLECTION_ID IN (", paste0("'", collection_IDs, "'", collapse = ","), ")") #Give me all fields from GEN_SAMPLED_FISH for the sillys in sillyvec
    
  dataAll0 <- RJDBC::dbGetQuery(con, qry)  #Pulling data from LOKI using the connection and fish data query
  
  discon <- RJDBC::dbDisconnect(con) # Disconnect from Loki
  
  dataAll <- dataAll0 %>% 
    dplyr::left_join(silly_collection_ID, by = "FK_COLLECTION_ID")# Adds silly code to ASL data
  
  if(import.vars == TRUE){
    
    output <- dataAll %>% 
      tibble::as_tibble() %>%
      dplyr::select(`Collection ID` = FK_COLLECTION_ID, 
                    `Fish ID` = FISH_ID, 
                    `Freshwater Age` = AGE_X, 
                    `Ocean Age` = AGE_Y, 
                    Sex = SEX, 
                    Length = LENGTH, 
                    Weight = WEIGHT, 
                    `Scale Card Number` = SCALE_CARD_NUM, 
                    `Scale Card Position` = SCALE_CARD_POS, 
                    `ASL Number` = ASL_NUMBER)
    
  } else{
    
    output <- dataAll %>% 
      tibble::as_tibble() %>%
      dplyr::select(`Silly code` = SILLY_CODE,
                           `Collection ID` = FK_COLLECTION_ID, 
                           `Fish ID` = FISH_ID, 
                           `Freshwater Age` = AGE_X, 
                           `Ocean Age` = AGE_Y, 
                           Sex = SEX, 
                           Length = LENGTH, 
                           Weight = WEIGHT, 
                           `Scale Card Number` = SCALE_CARD_NUM, 
                           `Scale Card Position` = SCALE_CARD_POS, 
                           `ASL Number` = ASL_NUMBER)
  }
  
  if(!is.null(file)){
    
    readr::write_csv(output, file, na = "", eol = "\r\n") # Write out a csv file without NAs
    
  }
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime)
  
  return(output)
  
}


#' @rdname get_asl_data
#' @export
ASL_Import.GCL <- get_asl_data  