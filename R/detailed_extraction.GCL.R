detailed_extraction.GCL <- function(plate_ids = NULL, sillyvec = NULL, username, password, file = NULL){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function connects to LOKI and pulls the well information for each DNA extraction plate similar to the Detailed Extraction tab in LOKI Collection Manager.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   plate_ids - a numeric vector of plate IDs (e.g. plate_ids <- c(4142, 7221, 7222, 8589))
  #
  #   sillyvec - optional, a character vector of silly codes (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART09")), used to get plate_ids if you do not know them
  #
  #   username - your state user name
  #
  #   password - your password used to access LOKI - see Eric Lardizabal if you don't have a password for LOKI
  #
  #   file - the file path, including .csv extension, for writing out a csv file of the output. 
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   The function outputs a tibble containing the following 5 variables:
  #  
  #   WELL_NO, FK_PLATE_ID, SILLY_CODE, FISH_NO, TISSUETYPE
  #
  #   If file is supplied, the tibble is written to a csv file with NAs removed.
  #   
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  #   password = "************"
  #
  #   source("~/../R/Functions.R")  # Load GCL functions; your path may differ.
  #
  #   detailed_extraction.GCL(plate_ids <- c(4142, 7221, 7222, 8589), sillyvec = NULL, username = "krshedd", password = password, file = "C:/Users/krshedd/Documents/R/test_wells_table.csv")
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #  This function was designed to pull all wells for each plate. 
  #  If you supply `sillyvec`, it will first query LOKI to figure out which plate_ids have those sillys, and then it will query LOKI again for all wells in those plate_ids.
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  
  if(!is.null(sillyvec) & !is.null(plate_ids)) {
    
    stop("You cannot supply both `sillyvec` and `plate_ids`, choose one or the other!")
    
  }  # error catching
  
  # This copies the "odbc8.jar" file to the R folder on your computer if it doesn't exist there. This file contains the java odbc drivers needed for RJDBC
  
  if(!file.exists(path.expand("~/R"))){
    
    dir <- path.expand("~/R")
    
    dir.create(dir)
    
    bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
    
  } else {
    
    if(!file.exists(path.expand("~/R/ojdbc8.jar"))){
      
      bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
      
    }
    
  }
  
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
    
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc8.jar", " ")  # https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
  } else {
    
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = path.expand("~/R/ojdbc8.jar"), " ")
    
  }
  
  url <- loki_url()  # This is a function that gets the correct URL to access the database on the oracle cloud
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)  # The database connection
  
  if(is.null(sillyvec)) {  # pulling by FK_PLATE_ID
    
    qry_plate_id <- paste("SELECT * FROM AKFINADM.V_GEN_DNA_WELL WHERE FK_PLATE_ID IN (", paste0("'", plate_ids, "'", collapse = ","), ")", sep = "")  # Query
    
  } else {  # pulling by SILLY_CODE, first get FK_PLATE_IDs, then pull by FK_PLATE_ID
    
    qry_sillyvec <- paste("SELECT * FROM AKFINADM.V_GEN_DNA_WELL WHERE SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")", sep = "")  # Query
    
    dataTmp <- RJDBC::dbGetQuery(con, qry_sillyvec)  # Pulling data from LOKI using the connection and sillyvec data query
    
    plate_ids <- unique(dataTmp$FK_PLATE_ID)  # Get all FK_PLATE_IDs
    
    qry_plate_id <- paste("SELECT * FROM AKFINADM.V_GEN_DNA_WELL WHERE FK_PLATE_ID IN (", paste0("'", plate_ids, "'", collapse = ","), ")", sep = "")  # Query
    
  }

  dataAll <- RJDBC::dbGetQuery(con, qry_plate_id)  # Pulling data from LOKI using the connection and wells table data query by FK_PLATE_ID
  
  discon <- RJDBC::dbDisconnect(con) # Disconnect from Loki
  
  output <- dataAll %>%
    tibble::as_tibble() %>%
    dplyr::select(
      WELL_NO,
      FK_PLATE_ID,
      SILLY_CODE,
      FISH_NO,
      TISSUETYPE
    )

  if(!is.null(file)){
    
    readr::write_csv(output, file, na = "") # Write out a csv file without NAs
    
  }
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime)
  
  return(output)}