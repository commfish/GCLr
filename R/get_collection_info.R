#' Get Collection Information
#' 
#' This function connects to LOKI and pulls a collection information report.
#' 
#' @param file The file path, including .csv extension, for writing out a csv file of the output. (default = NULL)
#' @param username Your state user name
#' @param password Your password used to access LOKI; see Eric Lardizabal if you need to set up a password
#' @param location.pmatch Logical; whether to use partial matching for location (default = FALSE)
#' @param sampled Pull only sampled (Yes) or unsampled (No) collection information. Choices: "Yes" or "No" (default = NULL); Variations of Yes and No will also work: "yes", "Y", "y", "no", "N", "n"
#' @param common_name The common name exactly the way it is spelled in Loki (e.g., "Salmon, Sockeye") (default = NULL)
#' @param collection_type Type of collection. Choices: "Baseline", "Experimental", "Forensics", "Mark Recapture", "Mixture", "N/A", "Species ID" (default = NULL)
#' @param region The region name exactly the way it is spelled in Loki (default = NULL)
#' @param quadrant The quadrant name exactly the way it is spelled in Loki (default = NULL)
#' @param location The location name exactly the way it is spelled in Loki (default = NULL)
#' 
#' @returns The function outputs a tibble of collection information and writes out an optional csv file
#' 
#' @note 
#' The counts in the 'Bulk remaining' column may not be accurate for some collections. 
#' This happens when tissues are transferred out of the bottles to another container type and the remaining bulk count is not updated in Loki.
#'       
#' @details
#' This function only requires username and password arguments to work. If a file path is supplied, then a csv file of the collection information will be written.
#' The remaining arguments, except location.pmatch, are filters. All filters can accept multiple filter items (i.e., vectors with lengths >= 1). If no filters are supplied then information for all collections in Loki will be pulled.
#' When supplying a region and quadrant, the quadrant must be contained within the region or the function will stop and list the allowable quadrant values for the region. 
#' When location.pmatch = FALSE, the location names must be spelled the same as in Loki. If a location name doesn't exist in Loki, the function will stop and an error message will list out the locations that don't exist.
#' When location.pmatch = TRUE, the function will find collections with location names containing text that matches the supplied location names. Using the location.pmatch argument can be useful if you don't know the exact spelling of the location names you are looking for.
#' If no collection location names contain the supplied location text string(s), the function will stop and give an error message.
#' 
#' @examples
#' \dontrun{
#' 
#' #All collection info and no output file
#' GCLr::get_collection_info(file = NULL, username = "awbarclay", password = "************", location.pmatch = FALSE, sampled = NULL, common_name = NULL, collection_type = NULL, region = NULL, quadrant = NULL, location = NULL)
#' 
#' #All collection info and write an output file
#' GCLr::get_collection_info(file = "CollectionInfo.csv", username = "awbarclay", password = "************", location.pmatch = FALSE, sampled = NULL, common_name = NULL, collection_type = NULL, region = NULL, quadrant = NULL, location = NULL)
#' 
#' #Collection info for all sockeye salmon baseline collections that were sampled
#' GCLr::get_collection_info(file = "CollectionInfo.csv", username = "awbarclay", password = "************", location.pmatch = FALSE, sampled = "Yes", common_name = "Salmon, Sockeye", collection_type = "Baseline", region = NULL, quadrant = NULL, location = NULL)
#' 
#' #Find all sampled Chinook salmon baseline collections with "Moose" in the collection location 
#' GCLr::get_collection_info(file = "CollectionInfo.csv", username = "awbarclay", password = "************", location.pmatch = TRUE, sampled = "Yes", common_name = "Salmon, Chinook", collection_type = "Baseline", region = NULL, quadrant = NULL, location = "Moose")
#' 
#' }
#' @export
get_collection_info <- function(file = NULL, username, password, location.pmatch = FALSE, sampled = NULL, common_name = NULL, collection_type = NULL, region = NULL, quadrant = NULL, location = NULL){
  
  options(java.parameters = "-Xmx4g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  # Database queries
  ## Collection table
  colqry <- paste("SELECT * FROM AKFINADM.GEN_COLLECTIONS") 
  
  col_All <- RJDBC::dbGetQuery(con, colqry) 
  
  ## Species table for common name
  spqry <- paste("SELECT * FROM AKFINADM.V_COLLECTIONS_SPECIES")#
  
  sp_All0 <- RJDBC::dbGetQuery(con, spqry) 
  
  ## Counts
  qry <- paste("SELECT * FROM AKFINADM.V_COLLECTION_COUNTS") 
  
  counts <- RJDBC::dbGetQuery(con, qry) 
  
  # Disconnect from database
  discon <- RJDBC::dbDisconnect(con) ## End database queries
  
  # Join collection and species data frames
  sp_All <- sp_All0 %>% 
    dplyr::select(SILLY_CODE, COMMON_NAME)
  
  col_info <- dplyr::left_join(col_All, sp_All, by = "SILLY_CODE") 
  
  # Join collection info with counts
  # Remove zeros from bulk and fish count columns before joining
  counts <- counts %>% 
    dplyr::select(COLLECTION_ID, BULK_REMAINING, FISH_COUNT) %>% 
    dplyr::mutate(BULK_REMAINING = dplyr::case_when(BULK_REMAINING == 0~NA,
                                             TRUE~BULK_REMAINING))%>% 
    dplyr::mutate(FISH_COUNT = dplyr::case_when(FISH_COUNT == 0~NA,
                                                    TRUE~FISH_COUNT))
  
  col_info <- dplyr::left_join(col_info, counts, by = "COLLECTION_ID")
  
  # Filter col_info data frame
  ## sampled
  if(!is.null(sampled)){
    
    if(sampled %in% c("y", "Y", "yes")){
      
      sampled <- "Yes"
      
    }
    
    if(sampled %in% c("n", "N", "no")){
      
      sampled <- "No"
      
    }
    
    if(!sampled %in% c("y", "Y", "yes", "Yes", "n", "N", "no", "No")){
      
      stop("Allowable values for the sampled arguement are 'Yes' or 'No'.")
      
    }
    
    col_info <- col_info %>% 
      dplyr::filter(SAMPLED == sampled)
    
  }
  
  ## common name
  if(!is.null(common_name)){
    
    if(any(!common_name %in% col_info$COMMON_NAME %>% unique())){
      
      stop("The supplied species common name does not match any common name in Loki. Please make sure the common name supplied is spelled and capitalized exactly the way it is in Loki.")
      
    }
    
    col_info <- col_info %>% 
      dplyr::filter(COMMON_NAME %in% common_name)
    
  }
  
  ## collection type
  if(!is.null(collection_type)){
    
    collection_type <- tools::toTitleCase(collection_type)# Make sure collection type is capitalized
    
    if(any(!collection_type %in% col_info$COLLECTION_TYPE %>% unique())){
      
      stop(paste0("Allowable values for collection type are ", paste0("'", col_info$COLLECTION_TYPE %>% unique(), collapse = "', "), "'"))
      
    }
    
    col_info <- col_info %>% 
      dplyr::filter(COLLECTION_TYPE %in% collection_type)
    
  }
  
  ## region
  if(!is.null(region)){
    
    region <- tools::toTitleCase(region)# Make sure region is capitalized
    
    if(any(!region %in% col_info$REGION_CODE %>% unique())){
      
      message <- paste0("Allowable values for region are: ", paste0("'", col_info$REGION_CODE %>% unique(), collapse = "', "), "'")
      
      options(warning.length = nchar(message)+1000)
      
      stop(message)
      
    }
    
    col_info <- col_info %>% 
      dplyr::filter(REGION_CODE %in% region)
    
  }
  
  ## quadrant
  if(!is.null(quadrant)){
    
    quadrant <- tools::toTitleCase(quadrant)# Make sure quadrant is capitalized
    
    if(any(!quadrant %in% col_info$QUADRANT %>% unique())){
      
      message <- paste0("Allowable values for quadrant are: ", paste0("'", col_info$QUADRANT %>% unique(), collapse = "', "), "'")
      
      options(warning.length = nchar(message)+1000)
        
      stop(message)
      
    }
    
    col_info <- col_info %>% 
      dplyr::filter(QUADRANT %in% quadrant)
    
  }
  
  ## location
  if(!is.null(location)){
   
   if(location.pmatch == TRUE){
     
     location0 <- sapply(location, function(loc){
       
       col_info$LOCATION_CODE[grepl(loc, col_info$LOCATION_CODE, ignore.case = TRUE)]
       
       }) %>% unlist() %>% unique()
     
     if(length(location0) == 0){
       
       stop("No collections had partial matching location names.")
       
     }else(location <- location0)
     
     col_info <- col_info %>% 
       dplyr::filter(LOCATION_CODE %in% location)
     
   } else{
     
     if(any(!location %in% col_info$LOCATION_CODE %>% unique())){
       
       nomatch <- setdiff(location, col_info$LOCATION_CODE)
       
       message <- paste0("The following location names are not in Loki:", paste0("'", nomatch, collapse = "', "), "'\nCheck your spelling or try setting loc.pmatch to TRUE")
       
       options(warning.length = nchar(message)+1000)
       
       stop(message)
       
       }
     
     col_info <- col_info %>% 
       dplyr::filter(LOCATION_CODE %in% location)
     
     }
    
    }else{
      
      if(location.pmatch == TRUE){
        
        warning("No location names supplied for partial matching.")
        
        }
      
      }## End filters
  

  # Clean up output
  output <- col_info %>% 
    dplyr::mutate(VERIFIED_DATE = lubridate::as_date(VERIFIED_DATE),
                  DATE_RECEIVED = lubridate::as_date(DATE_RECEIVED),
                  COLLECTION_DATE = lubridate::as_date(COLLECTION_DATE)) %>% 
    dplyr::select(`Collection Type` = COLLECTION_TYPE, 
           Quality = QUALITY, 
           `Life Stage` = LIFE_STAGE, 
           `Common Name` = COMMON_NAME,
           `Collection ID` = COLLECTION_ID, 
           `Silly Code` = SILLY_CODE, 
           `Alternate ID` = COL_ID, 
           `Sampled` = SAMPLED, 
           `Bulk remaining` = BULK_REMAINING,
           `Fish count` = FISH_COUNT,
           Status = STATUS, 
           `Satus Comment` = STATUS_COMMENT, 
           `Permit Number` = PERMIT_NUMBER, 
           `AWC Stream Code` = AWC_STREAM_CODE, 
           `Location Descriptor` = LOCATION_DESCRIPTOR, 
           Region = REGION_CODE, 
           Quadrant = QUADRANT, 
           Location = LOCATION_CODE, 
           Latitude = LATTITUDE, 
           Longitude = LONGITUDE, 
           `Collection Date` = COLLECTION_DATE, 
           `Date Received` = DATE_RECEIVED, 
           `Is Verified` = VERIFIED_YN, 
           `Verified By` = VERIFIED_BY, 
           `Verified Date` = VERIFIED_DATE, 
           `Collection Comment` = COL_COMMENT)
  
  
  # If file supplied, write out csv
  if(!is.null(file)){
    
    readr::write_csv(output, file = file, na = "")
    
  }
  
  return(output)

}