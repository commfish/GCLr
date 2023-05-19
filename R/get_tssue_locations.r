get_tssue_locations <- function(unit, username, password, bad_locations = TRUE, all_tissues = TRUE) {
  #########################################
  # This function pulls tissue data from OceanAK, for a given storage location, and updates the tissue location maps.
  # Created by: Chase Jalbert
  # Created on: 9/1/2020
  #
  # Inputs~~~~~~~~~~~~~~~~~~
  #  unit - the storage unit(s) you're interested in; can be single or multiple units (e.g., "99" or c("WA","WB","WC"))
  #         Note - this works best with a single type of unit (e.g., all of B7, all of Warehouse). You can combine but the output could get messy due to different shelf naming conventions
  #  username - your LOKI username and password
  #  password - your LOKI password
  #  bad_locations - TRUE or FALSE, do you want a list of tissues with incorrect location information; default = FALSE
  #  all_tissues- TRUE or FALSE, do you want a list of ALL tissues [potentially large file]]; default = FALSE
  #
  #
  # Outputs~~~~~~~~~~~~~~~~~
  #  tissuemap - an object mapping out all the tissues for a storage unit
  #  all_tissues- an object consisting of all tissues within OceanAK, within a storage unit
  #  bad_locations - an object consisting of all tissues within OceanAK, containing unexpected, missing, or otherwise incorrect location information (i.e., SHELF_RACK and/or SLOT)
  #
  #
  # Example~~~~~~~~~~~~~~~~~
  #  username = "awesomeuser"
  #  password = "awesomepassword1"
  #  unit = c("A", "B", "C", "D", "E", "F") # We want these units in B7
  #  bad_locations = TRUE # yes, I want to identify all wrong/missing locations
  #  all_tissues = TRUE # yes, I want all the data
  #
  #  TissueLocations_2R.GCL(unit = unit, username = username, password = password, bad_location = bad_location, all_data = all_data)
  #  
  # 
  ##########################################
  
  # Setup  
  
  # Source R functions:
  source("C:/Users/hahoyt1/Documents/R/Functions.GCL.R")
  
  # Database setup
  ## This copies the "odbc8.jar" file to the R folder on your computer if it doesn't exist there. This file contains the java odbc drivers needed for RJDBC
  if (!file.exists(path.expand("~/R"))) {
    dir <- path.expand("~/R")
    
    dir.create(dir)
    
    bool <-
      file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
    
  } else {
    if (!file.exists(path.expand("~/R/ojdbc8.jar"))) {
      bool <-
        file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
      
    }
    
  }
  
  start.time <- Sys.time()
  
  options(java.parameters = "-Xmx10g")
  
  if (file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
    drv <-
      RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc8.jar", " ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
  } else {
    drv <-
      RJDBC::JDBC("oracle.jdbc.OracleDriver",
                  classPath = path.expand("~/R/ojdbc8.jar"),
                  " ")
    
  }
  
  ## Build database URL
  url <-
    loki_url() # This is a function that gets the correct URL to access the database on the oracle cloud
  
  ## Connect to database
  con <-
    RJDBC::dbConnect(
      drv = drv,
      url = url,
      user = username,
      password = password
    ) #The database connection
  
  ## Create the query
  gnoqry <-
    paste(
      "SELECT * FROM AKFINADM.V_GEN_SAMPLED_FISH_TISSUE WHERE UNIT IN (",
      paste0("'", unit, "'", collapse = ","),
      ")",
      sep = ""
    )
  
  # Data import
  ## Open the connection and pull data from the database
  dataAll <-
    RJDBC::dbGetQuery(conn = con, statement = gnoqry)
  
  ## Subset the data
  dataSubset <- dataAll %>%
    dplyr::filter(
      !grepl("ERICTEST", SILLY_CODE, ignore.case = TRUE),
      # drop Eric's test sillys
      PK_TISSUE_TYPE != "DNA"
      # drop DNA tissue type
      #is.na(EXHAUSTED_HOW) |
      #  EXHAUSTED_HOW != "Discarded"
    ) %>% # drop tissues marked as discarded
    dplyr::select(
      c(
        FK_COLLECTION_ID,
        SILLY_CODE,
        FK_FISH_ID,
        STORAGE_ID,
        PK_TISSUE_TYPE,
        UNIT,
        SHELF_RACK,
        SLOT,
        EXHAUSTED_HOW
      )
    ) %>% # select columns
    tidyr::unite(
      col = "tissue_id",
      SILLY_CODE,
      PK_TISSUE_TYPE,
      STORAGE_ID,
      sep = "_",
      remove = FALSE
    ) # create unique tissue identifier
  
  
  # Identify incorrect or missing tissue locations
  ## Note - SHELF_RACK and SLOT terminology differs with each storage location, so I've broken them out below
  bad_location <- dataSubset %>%
    tidyr::unite(col = "shelf_id", SHELF_RACK, SLOT, remove = FALSE) %>% # combine these into single column for an id
    # Find anything that does NOT match acceptable storage location code, for each unit
    dplyr::mutate(
      wrong = dplyr::case_when(
        # Freezer 99 - acceptable locations are: ###_UppercaseLetter
        UNIT %in% unit[stringr::str_detect(unit, "99")] &
          !stringr::str_detect(string = shelf_id, pattern = "[0-9]{3}_[A-Z]") ~ "wrong",
        # Warehouse - acceptable locations are: W[1-33]_[1-10]
        UNIT %in% unit[stringr::str_detect(unit, "^(W[A-Z])$")] &
          !str_detect(string = shelf_id, pattern = "^([1-9]|1[0-9]|2[0-9]|3[0-3])_([1-9]|10)$") ~ "wrong",
        # Freezer 7 or 8 - acceptable locations are: UppercaseLetter[A-E]_[1-6]
        UNIT %in% unit[stringr::str_detect(unit, "^([7-8])$")] &
          !stringr::str_detect(string = shelf_id, pattern = "^([A-E])_([1-6])$") ~ "wrong",
        # B7 - acceptable locations are:  [1-10]_UppercaseLetter
        UNIT %in% unit[stringr::str_detect(unit, "^([A-W])$")] &
          !stringr::str_detect(string = shelf_id, pattern = "^([1-9]|10)_([A-Z])$") ~ "wrong",
        # Freezer 1 - acceptable locations are: UppercaseLetter[A-E]_[1-6]
        UNIT %in% unit[stringr::str_detect(unit, "1")] &
          !stringr::str_detect(string = shelf_id, pattern = "^([A-E]_[1-6])$") ~ "wrong",
        TRUE ~ "correct"
      )
    ) %>%
    dplyr::filter(wrong == "wrong") %>%
    dplyr::select(-wrong)
  
  ## If, you said TRUE in setup, then assign the bad location object to your environment. From here you can export or do whatever you want with it.
  if (bad_locations == TRUE) {
    
    # export CSV of incorrect tissues
    write_csv(x = bad_location, path = paste0("V:/Lab/Archive Storage/Archive Sample Maps from R/bad_tissue_locations_", format(Sys.Date(), format = "%Y%m%d"), ".csv")) 
    
    assign(
      x = "bad_location",
      value = bad_location,
      pos = 1,
      envir = .GlobalEnv
    )
    message(paste(
      "All tissues with incorrect locations stored in object `bad_location` and a CSV was output"))
  }
  
  # Create final map
  ## Make fish number ranges
  fish_range <- dataSubset %>%
    dplyr::anti_join(
      bad_location,
      by = c(
        "FK_COLLECTION_ID",
        "tissue_id",
        "SILLY_CODE",
        "FK_FISH_ID",
        "STORAGE_ID",
        "PK_TISSUE_TYPE",
        "UNIT",
        "SHELF_RACK",
        "SLOT"
      )
    ) %>% # dropping any unknown/incorrect location collections
    dplyr::group_by(tissue_id) %>% # for each tissue id:
    dplyr::summarise(min = min(FK_FISH_ID), # find the min fish number
                     max = max(FK_FISH_ID), # find the max fish number
                     .groups = "drop_last") 
  
  ## Assign the fish range to each fish
  TissueData <-
    fish_range %>% # contains tissue range and does not include bad_location tissues
    dplyr::left_join(
      dataSubset %>% # Original data
        dplyr::anti_join(
          bad_location,
          by = c(
            "FK_COLLECTION_ID",
            "tissue_id",
            "SILLY_CODE",
            "FK_FISH_ID",
            "STORAGE_ID",
            "PK_TISSUE_TYPE",
            "UNIT",
            "SHELF_RACK",
            "SLOT"
          )
        ) %>% # dropping any unknown/incorrect location collections
        dplyr::select(tissue_id, UNIT, SHELF_RACK, SLOT) %>% # just grab a few columns from this data set
        tidyr::unite(col = "shelf_id", SHELF_RACK, SLOT, remove = FALSE),
      # make shelf location info
      by = "tissue_id"
    ) %>% # join fish_range with dataSubset by tissue_id, keeping only entries in fish_range
    dplyr::distinct() # grab just unique values
  
  
  ## Convert into map format
  data_map <- TissueData %>%
    tidyr::unite(range, min, max, sep = "-") %>%
    tidyr::separate(tissue_id,
                    into = c("silly", "tissue", "barcode"),
                    sep = "_") %>% # we have to split these to insert fish range
    tidyr::unite(map_id, silly, tissue, range, barcode, sep = " ") %>% # now making the final map entries
    dplyr::distinct() %>% # we just want unique entries
    dplyr::select(-c(shelf_id)) %>% # drop extra column
    dplyr::arrange(UNIT, SHELF_RACK, SLOT) %>% # sorting by shelf_rack, then slot
    dplyr::group_by(UNIT, SHELF_RACK, SLOT) %>%
    dplyr::group_modify( ~ {
      .x %>%
        dplyr::mutate(row = row_number())
    }) # experimental magic... - modifies groups, based on previous group_by(); I am simply getting row numbers here
  
  ## Create tissue location map
  tissuemap <- data_map %>%
    tidyr::pivot_wider(names_from = SLOT,
                       values_from = map_id) %>% # move from long to wide format
    dplyr::select(-row) %>% # ditch the row ids from earlier
    dplyr::select(UNIT, SHELF_RACK, sort(tidyselect::peek_vars())) # reorder into UNIT, SHELF_RACK, anything else (i.e., shelves order)

  
  ## output tissue map to environment for saving, always do this
  assign(
    x = "tissuemap",
    value = tissuemap,
    pos = 1,
    envir = .GlobalEnv
  )
  
  ## write a csv of the tissue map to the folder:
  write_csv(x = tissuemap, path = paste0("V:/Lab/Archive Storage/Archive Sample Maps from R/tissuemap_", format(Sys.Date(), format = "%Y%m%d"), ".csv"))
  
 # Wrap up function 
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  
  print(fulltime)
  message(paste0("Map of tissue locations stored in object 'tissuemap' and a CSV was output"))
  
  ## If, you said TRUE in setup, then assign the all_data object to your environment. From here you can export or do whatever you want with it.
  if (all_tissues == TRUE) {
    
    # export CSV of ALL tissues
    write_csv(x = dataAll, path = paste0("V:/Lab/Archive Storage/Archive Sample Maps from R/all_tissues_", format(Sys.Date(), format = "%Y%m%d"), ".csv"))
    
    assign(
      x = "all_tissues",
      value = dataAll,
      pos = 1,
      envir = .GlobalEnv
    )
    
    message(paste("All fish (raw) from database stored in object `all_data` and a CSV was output"))
  }
}