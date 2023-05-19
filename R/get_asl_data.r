get_asl_data <- function(sillyvec, username, password, file = NULL, import.vars = TRUE){
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function connects to Loki and pulls the fish information (ASL data) for each silly in sillyvec.
  #  This pulls the same information as the "ASL Import" report in OceanAK.
  #  
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   sillyvec - a character vector of silly codes
  #
  #   username - your state user name
  #
  #   password - your password used to access Loki - see Eric Lardizabal if you don't have a password for LOKI
  #
  #   file - the file path, including .csv extension, for writing out a csv file of the output. 
  #
  #   import.vars - logical; do you only want the variables allowed for the ASL importer?
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   If import.vars = TRUE, the function outputs a tibble containing the following variables: 
  #   Collection ID, Fish ID, Freshwater Age, Ocean Age, Sex, Length, Weight, Scale Card Number, Scale Card Position, ASL Number
  #
  #   If import.vars = FALSE, a variable for silly code is added to the output tibble.
  #
  #   If file is supplied, the tibble is written to a csv file with NAs removed.
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  #   password = "************"
  #
  #   source("~/R/Functions.R")# Load GCL functions; your path may differ.
  #
  #   get_asl_data(sillyvec = c("PPORTCH21", "PBARABCR21"), username = "awbarclay", password = password, import.vars = FALSE, file = "C:/Users/awbarclay/Documents/R/test_fish_table.csv") #This example includes silly code in output.
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # The output of this function is can be used for creating an import file for the Loki ASL Data Importer. 
  #
  # The the importer requires that all columns are in the correct order and spelled correctly. The importer will not work if you include silly code in the import file.
  #
  # When writing the import file, make sure it does not contain NAs. (see example below) 
  #
  # Also, when using readr::write_csv to write out the tissue import file, make sure to change the eol argument to from the default "\n" to "\r\n" or the importer
  # will give you an error message about the header names.  e.g.,  import_file %>% write_csv(file = "ImportFile.csv", na = "", eol = "\r\n")
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  
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
    
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc8.jar", " ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
  } else {
    
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = path.expand("~/R/ojdbc8.jar"), " ")
    
  }
  
  url <- loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password) #The database connection
  
  qry <- paste("SELECT * FROM AKFINADM.V_GEN_SAMPLED_FISH_TISSUE WHERE SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")", sep = "") #Query the tissues table to get collection IDs. Need to see if Eric can add silly code to the fish table.
  
  silly_collection_ID <- RJDBC::dbGetQuery(con, qry) %>% 
    select(SILLY_CODE, FK_COLLECTION_ID) %>% 
    distinct()#This gets the collection IDs.
  
  collection_IDs <- silly_collection_ID$FK_COLLECTION_ID
  
  qry <- paste0("SELECT * FROM AKFINADM.GEN_SAMPLED_FISH WHERE FK_COLLECTION_ID IN (", paste0("'", collection_IDs, "'", collapse = ","), ")") #Give me all fields from GEN_SAMPLED_FISH for the sillys in sillyvec
    
  dataAll0 <- RJDBC::dbGetQuery(con, qry)  #Pulling data from LOKI using the connection and fish data query
  
  discon <- RJDBC::dbDisconnect(con) # Disconnect from Loki
  
  dataAll <- dataAll0 %>% 
    left_join(silly_collection_ID, by = "FK_COLLECTION_ID")# Adds silly code to ASL data
  
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
    
    write_csv(output, file, na = "", eol = "\r\n") # Write out a csv file without NAs
    
  }
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime)
  
  return(output)
  
}