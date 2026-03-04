#' Get DNA Extraction Information
#' 
#' This function connects to LOKI and pulls the well information for each DNA extraction plate.
#' 
#' @param plate_ids a numeric vector of plate IDs (e.g. plate_ids <- c(4142, 7221, 7222, 8589)) (default = NULL)
#' @param sillyvec optional, a character vector of silly codes (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART09")), used to get plate_ids if you do not know them (default = NULL)
#' @param password - your password used to access LOKI - contact Eric Lardizabal if you don't have a password for LOKI
#' @param file the file path, including .csv extension, for writing out a csv file of the output. (default = NULL)
#' 
#' @returns The function outputs a tibble containing the following 8 variables:
#'  \itemize{
#'         \item \code{WELL_NO}: 
#'         \item \code{FK_PLATE_ID}: the extraction plate ID
#'         \item \code{SILLY_CODE}: the collection silly code
#'         \item \code{FISH_NO}: the collection fish number
#'         \item \code{TISSUETYPE}: e.g., Axillary Process, Fin, Heart, etc.
#'         \item \code{FREEZER}: The freezer number
#'         \item \code{RACK}: The freezer rack
#'         \item \code{SLOT}: The freezer slot
#'         
#'       }
#'       
#' @details
#' This function can pull extraction information from Loki using either the plate IDs (if known) or using silly codes (if plate IDs aren't known).  
#' This function was designed to pull all wells for each plate. If you supply `sillyvec`, it will first query LOKI to figure out which plate_ids have those sillys, and then it will query LOKI again for all wells in those plate_ids.
#' When using plate IDs, the output will only contain extraction information for the supplied plate IDs. 
#' If using silly codes, the output will contain information for all extraction plates for the supplied silly codes. 
#' 
#' @examples
#' \dontrun{
#'  
#' password = "************"
#' 
#' #By plate ID
#' GCLr::get_extraction_info(plate_ids = c(4142, 7221, 7222, 8589), sillyvec = NULL, username = "krshedd", password = password, file = path.expand("~/plateID_extraction_info.csv") )
#' 
#' #By silly code
#' GCLr::get_extraction_info(plate_ids = NULL, sillyvec = c("KKENM06", "KQUART06", "KCRESC06", "KRUSSR08", "KQUART08", "KDAVE08", "KSLIK08", "KKENU09", "KQUART09", "KKENN09"), username = "krshedd", password = password, file = path.expand("~/silly_extraction_info.csv") )
#' 
#' }
#' @export
get_extraction_info <- function(plate_ids = NULL, sillyvec = NULL, username, password, file = NULL){
  
  if(!is.null(sillyvec) & !is.null(plate_ids)) {
    
    stop("You cannot supply both `sillyvec` and `plate_ids`, choose one or the other!")
    
  }  # error catching
  
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  if(is.null(sillyvec)) {  # pulling by FK_PLATE_ID
    
    qry_plate_id_fish <- paste("SELECT * FROM AKFINADM.V_GEN_DNA_WELL WHERE FK_PLATE_ID IN (", paste0("'", plate_ids, "'", collapse = ","), ")", sep = "")  # Query
    
  } else {  # pulling by SILLY_CODE, first get FK_PLATE_IDs, then pull by FK_PLATE_ID
    
    qry_sillyvec <- paste("SELECT * FROM AKFINADM.V_GEN_DNA_WELL WHERE SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")", sep = "")  # Query
    
    dataTmp <- RJDBC::dbGetQuery(con, qry_sillyvec)  # Pulling data from LOKI using the connection and sillyvec data query
    
    plate_ids <- unique(dataTmp$FK_PLATE_ID)  # Get all FK_PLATE_IDs
    
    qry_plate_id_fish <- paste("SELECT * FROM AKFINADM.V_GEN_DNA_WELL WHERE FK_PLATE_ID IN (", paste0("'", plate_ids, "'", collapse = ","), ")", sep = "")  # Query
    
  }
  
  # Query for extraction plate locations
  qry_plate_id_plate <- paste("SELECT * FROM AKFINADM.GEN_DNA_EXTRACTION WHERE PLATE_ID IN (", paste0("'", plate_ids, "'", collapse = ","), ")", sep = "")  # Query the table containing freezer locations for each plate ID
  
  
  extraction_info_by_individual0 <- RJDBC::dbGetQuery(con, qry_plate_id_fish)  # Pulling extraction data by fish ID
  
  extraction_info_by_plate0 <- RJDBC::dbGetQuery(con, qry_plate_id_plate) # Pulling extraction data by plate ID

  
  discon <- RJDBC::dbDisconnect(con) # Disconnect from Loki
  
  extraction_info_by_individual <- extraction_info_by_individual0 %>%
    tibble::as_tibble() %>%
    dplyr::select(
      WELL_NO,
      FK_PLATE_ID,
      SILLY_CODE,
      FISH_NO,
      TISSUETYPE
    )
  
  extraction_info_by_plate <- extraction_info_by_plate0 %>% 
    tibble::as_tibble() %>%
    dplyr::select(FK_PLATE_ID = PLATE_ID, FREEZER, RACK, SLOT)
  
  output <- left_join(extraction_info_by_individual, extraction_info_by_plate, by = "FK_PLATE_ID")
  
  if(!is.null(file)){
    
    readr::write_csv(output, file, na = "") # Write out a csv file without NAs
    
  }
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime)
  
  return(output)
  
}