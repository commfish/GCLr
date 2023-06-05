get_tissue_data <- function(sillyvec, username, password, file = NULL, import.vars = TRUE){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function connects to LOKI and pulls the tissue information for each silly in sillyvec.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   sillyvec - a character vector of silly codes (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART09"))
  #
  #   username - your state user name
  #
  #   password - your password used to access LOKI - see Eric Lardizabal if you don't have a password for LOKI
  #
  #   file - the file path, including .csv extension, for writing out a csv file of the output. 
  #
  #   import.vars - if TRUE (default), the output will only contain the 31 fields used by the Loki tissue importer (same as OceanAK report)
  #                 if FALSE, the output will include 7 additional collection information fields
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   If import.vars = TRUE, the function outputs a tibble containing the following 31 variables:
  #  
  #   FK_COLLECTION_ID, FK_FISH_ID, PK_TISSUE_TYPE, CAPTURE_LOCATION, CAPTURE_DATE, END_CAPTURE_DATE, LATITUDE, LONGITUDE, MESH_SIZE, 
  #   MESH_SIZE_COMMENT, IS_MISSING_PAIRED_DATA_EXISTS, WELL_HAS_MORE_THAN_ONE_SAMPLE, IS_PRESENT_IN_DATASHEET, IS_PRESENT_BUT_NOT_IN_DS, 
  #   VIAL_BARCODE, CONTAINER_ARRAY_TYPE_ID, DNA_TRAY_WORKBENCH_ID, DNA_TRAY_CODE, DNA_TRAY_WELL_POS, DNA_TRAY_WELL_CODE, STORAGE_ID, 
  #   UNIT, SHELF_RACK, SLOT, EXHAUSTED_HOW, EXHAUSTED_BY, EXHAUSTED_DATE, AGENCY, OTHER_AGENCY_KEY, NUM_OTOLITHS_MISSING, OTO_INVENTORY_COMMENT
  #
  #   If import.vars = FALSE, the function outputs a tibble containing 7 collection information variables (columns 1-7) and the 31 tissue import variables (columns 8-38). 
  #   Here are the collection information variables:
  #
  #   SILLY_CODE, REGION_CODE, QUADRANT, LOCATION_CODE, LOCATION_DESCRIPTOR, LIFE_STAGE, COLLECTION_TYPE
  #
  #   If file is supplied, the tibble is written to a csv file with NAs removed.
  #   
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  #   password = "************"
  #
  #   source("~/R/Functions.R")# Load GCL functions; your path may differ.
  #
  #   get_tissue_data(sillyvec = c("KCDVF18", "KCDVF19", "KCDVF20"), username = "awbarclay", password = password, file = "C:/Users/awbarclay/Documents/R/test_tissue_table.csv", import.vars = FALSE)
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # The output of this function is can be used for creating an import file for the Loki tissue importer to add additional information to Loki tissue table (default output) and 
  # and extraction lists for lab projects.
  #
  # The tissue importer requires that all columns are in the correct order and spelled correctly. 
  #
  # When writing the import file, make sure it does not contain NAs. (see example above) 
  #
  # The tissue importer doesn't like some date formats. Before importing tissue information with dates, 
  # open the file in excel and change the format of date columns to numeric and then save the file; the numbers will get converted to dates in Loki.
  #
  # Also, when using readr::write_csv to write out the tissue import file, make sure to change the eol argument to from the default "\n" to "\r\n" or the importer
  # will give you an error message about the header names.  e.g.,  import_file %>% write_csv(file = "ImportFile.csv", na = "", eol = "\r\n")
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
 
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  qry <- paste("SELECT * FROM AKFINADM.V_GEN_SAMPLED_FISH_TISSUE WHERE SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")", sep = "") #Query
  
  dataAll <- RJDBC::dbGetQuery(con, qry)  #Pulling data from LOKI using the connection and tissue data query
  
  discon <- RJDBC::dbDisconnect(con) # Disconnect from Loki
  
  if(import.vars == TRUE){
    
    output <- dataAll %>% 
      tibble::as_tibble() %>%
      dplyr::mutate(CAPTURE_DATE = lubridate::as_date(CAPTURE_DATE), END_CAPTURE_DATE = lubridate::as_date(END_CAPTURE_DATE)) %>% 
      dplyr::select(FK_COLLECTION_ID, FK_FISH_ID, PK_TISSUE_TYPE, CAPTURE_LOCATION, CAPTURE_DATE, END_CAPTURE_DATE, LATITUDE, LONGITUDE, MESH_SIZE, 
             MESH_SIZE_COMMENT, IS_MISSING_PAIRED_DATA_EXISTS, WELL_HAS_MORE_THAN_ONE_SAMPLE, IS_PRESENT_IN_DATASHEET, IS_PRESENT_BUT_NOT_IN_DS, 
             VIAL_BARCODE, CONTAINER_ARRAY_TYPE_ID, DNA_TRAY_WORKBENCH_ID, DNA_TRAY_CODE, DNA_TRAY_WELL_POS, DNA_TRAY_WELL_CODE, STORAGE_ID, 
             UNIT, SHELF_RACK, SLOT, EXHAUSTED_HOW, EXHAUSTED_BY, EXHAUSTED_DATE, AGENCY, OTHER_AGENCY_KEY, NUM_OTOLITHS_MISSING, OTO_INVENTORY_COMMENT)
    
  } else{
    
    output <- dataAll %>% 
      tibble::as_tibble() %>%
      dplyr::mutate(CAPTURE_DATE = lubridate::as_date(CAPTURE_DATE), END_CAPTURE_DATE = lubridate::as_date(END_CAPTURE_DATE)) %>% 
      dplyr::select(SILLY_CODE, REGION_CODE, QUADRANT, LOCATION_CODE, LOCATION_DESCRIPTOR, LIFE_STAGE, COLLECTION_TYPE,
             FK_COLLECTION_ID, FK_FISH_ID, PK_TISSUE_TYPE, CAPTURE_LOCATION, CAPTURE_DATE, END_CAPTURE_DATE, LATITUDE, LONGITUDE, MESH_SIZE, 
             MESH_SIZE_COMMENT, IS_MISSING_PAIRED_DATA_EXISTS, WELL_HAS_MORE_THAN_ONE_SAMPLE, IS_PRESENT_IN_DATASHEET, IS_PRESENT_BUT_NOT_IN_DS, 
             VIAL_BARCODE, CONTAINER_ARRAY_TYPE_ID, DNA_TRAY_WORKBENCH_ID, DNA_TRAY_CODE, DNA_TRAY_WELL_POS, DNA_TRAY_WELL_CODE, STORAGE_ID, 
             UNIT, SHELF_RACK, SLOT, EXHAUSTED_HOW, EXHAUSTED_BY, EXHAUSTED_DATE, AGENCY, OTHER_AGENCY_KEY, NUM_OTOLITHS_MISSING, OTO_INVENTORY_COMMENT)
    
  }
  
  if(!is.null(file)){
    
    write_csv(output, file, na = "") # Write out a csv file without NAs
    
  }
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime)
  
  return(output)
  
}