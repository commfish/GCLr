#' Get Tissue Data from Loki
#' 
#' This function pulls tissue data from Loki, for a given storage location, and updates the tissue location maps.
#'
#' @param unit the storage unit(s) you're interested in; can be single or multiple units (e.g., "99" or c("WA","WB","WC")) (see details)
#'         
#' @param username your Loki username
#' @param password your Loki password
#' @param bad_locations TRUE or FALSE, do you want a list of tissues with incorrect location information (default = FALSE)
#' @param all_tissues  TRUE or FALSE, do you want a list of ALL tissues (default = FALSE) (see details)
#' 
#' @returns this function assigns the following objects to your workspace:
#' \itemize{
#'    \item \code{tissuemap}: an object mapping out all the tissues for a storage unit
#'    \item \code{all_tissues}: an object consisting of all tissues Loki, within a storage unit
#'    \item \code{bad_locations}: an object consisting of all tissues within Loki, containing unexpected, missing, or otherwise incorrect location information (i.e., SHELF_RACK and/or SLOT)
#'    }
#' @details
#' This function best when run with a single type of unit (e.g., all of B7, all of Warehouse). 
#' You can combine but the output could get messy due to different shelf naming conventions.
#' If you set `all_tissues = TRUE` the resulting output will be very large!
#'    
#'    
#' @examples
#' \dontrun{
#'  unit <- c("A", "B", "C", "D", "E", "F") # We want these units in B7
#'  
#'  get_tissue_locations(unit = unit, username = "awbarclay", password = password, bad_locations = TRUE, all_data = TRUE)
#' 
#' }
#' 
#' @export
get_tissue_locations <- function(unit, username, password, bad_locations = TRUE, all_tissues = TRUE) {
  
  start.time <- Sys.time()
  
  # Create database connection
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  ## Create the query
  gnoqry <-paste("SELECT * FROM AKFINADM.V_GEN_SAMPLED_FISH_TISSUE WHERE UNIT IN (", paste0("'", unit, "'", collapse = ","), ")", sep = "")
  
  # Data import
  ## Open the connection and pull data from the database
  dataAll <- RJDBC::dbGetQuery(conn = con, statement = gnoqry)
  
  discon <- RJDBC::dbDisconnect(con) #Disconnect from database
  
  ## Subset the data
  dataSubset <- dataAll %>% 
    dplyr::filter(!grepl("ERICTEST", SILLY_CODE, ignore.case = TRUE), # drop Eric's test sillys
                  PK_TISSUE_TYPE != "DNA") %>% # drop tissues marked as discarded. These were in the filter but commented out: is.na(EXHAUSTED_HOW) | EXHAUSTED_HOW != "Discarded"
    dplyr::select(
        FK_COLLECTION_ID,
        SILLY_CODE,
        FK_FISH_ID,
        STORAGE_ID,
        PK_TISSUE_TYPE,
        UNIT,
        SHELF_RACK,
        SLOT,
        EXHAUSTED_HOW
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
          !stringr::str_detect(string = shelf_id, pattern = "^([1-9]|1[0-9]|2[0-9]|3[0-3])_([1-9]|10)$") ~ "wrong",
        # Freezer 101,102,103 - acceptable locations are: UppercaseLetter[A-E]_[1-6]
        UNIT %in% unit[stringr::str_detect(unit, "101|102|103")] &
          !stringr::str_detect(string = shelf_id, pattern = "^([A-E])_([1-6])$") ~ "wrong",
        # B7 - acceptable locations are:  [1-10]_UppercaseLetter
        UNIT %in% unit[stringr::str_detect(unit, "^([A-W])$")] &
          !stringr::str_detect(string = shelf_id, pattern = "^([1-9]|10)_([A-Z])$") ~ "wrong",
        TRUE ~ "correct"
      )
    ) %>%
    dplyr::filter(wrong == "wrong") %>%
    dplyr::select(-wrong)
  
  ## If, you said TRUE in setup, then assign the bad location object to your environment. From here you can export or do whatever you want with it.
  if (bad_locations == TRUE) {
    
    # export CSV of incorrect tissues
    readr::write_csv(x = bad_location, file = paste0("V:/Lab/Archive Storage/Archive Sample Maps from R/bad_tissue_locations_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")) 
    
    assign(
      x = "bad_location",
      value = bad_location,
      pos = 1,
      envir = .GlobalEnv
    )
    message(paste(
      "All tissues with incorrect locations stored in object `bad_location` and a CSV was output to V:/Lab/Archive Storage/Archive Sample Maps from R/"))
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
    tidyr::unite(map_id, barcode, silly, tissue, range,  sep = " ") %>% # now making the final map entries
    dplyr::distinct() %>% # we just want unique entries
    dplyr::select(-c(shelf_id)) %>% # drop extra column
    dplyr::arrange(UNIT, SHELF_RACK, SLOT, map_id) %>% # sorting by shelf_rack, then slot
    dplyr::group_by(UNIT, SHELF_RACK, SLOT) %>%
    dplyr::group_modify( ~ {
      .x %>%
        dplyr::mutate(row = dplyr::row_number())
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
  readr::write_csv(x = tissuemap, file = paste0("V:/Lab/Archive Storage/Archive Sample Maps from R/tissuemap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
  
 # Wrap up function 
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime)
  message(paste0("Map of tissue locations stored in object 'tissuemap' and a CSV was output to V:/Lab/Archive Storage/Archive Sample Maps from R/"))
  
  ## If, you said TRUE in setup, then assign the all_data object to your environment. From here you can export or do whatever you want with it.
  if (all_tissues == TRUE) {
    
    # export CSV of ALL tissues
    readr::write_csv(x = dataAll, file = paste0("V:/Lab/Archive Storage/Archive Sample Maps from R/all_tissues_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
    
    assign(
      x = "all_tissues",
      value = dataAll,
      pos = 1,
      envir = .GlobalEnv
    )
    
    message(paste("All fish (raw) from database stored in object `all_data` and a CSV was output to V:/Lab/Archive Storage/Archive Sample Maps from R/"))
  
    }
}
