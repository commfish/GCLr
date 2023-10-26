#' Get Tissue Data
#'
#' This function connects to LOKI and pulls the tissue information for each silly in sillyvec.
#'
#' @param sillyvec A character vector of silly codes
#' @param unit the archive storage unit
#' @param shelf.rack the archive shelf.rack
#' @param username Your state user name.
#' @param password Your password used to access LOKI.
#' @param file The file path for writing out a CSV file of the output (optional).
#' @param import.vars If TRUE, output contains the 31 fields used by the Loki tissue importer (same as OceanAK report).
#' If FALSE, output includes 7 additional collection information fields, including SILLY_CODE.
#'
#' @returns A tibble containing the specific tissue data, as determined by `import.vars`.
#' 
#' @section Returns:
#' If \code{import.vars = TRUE}, the function returns a tibble with the following 31 variables:
#' \itemize{
#'   \item{FK_COLLECTION_ID}{...}
#'   \item{FK_FISH_ID}{...}
#'   \item{PK_TISSUE_TYPE}{...}
#'   \item{CAPTURE_LOCATION}{...}
#'   \item{CAPTURE_DATE}{...}
#'   \item{END_CAPTURE_DATE}{...}
#'   \item{LATITUDE}{...}
#'   \item{LONGITUDE}{...}
#'   \item{MESH_SIZE}{...}
#'   \item{MESH_SIZE_COMMENT}{...}
#'   \item{IS_MISSING_PAIRED_DATA_EXISTS}{...}
#'   \item{WELL_HAS_MORE_THAN_ONE_SAMPLE}{...}
#'   \item{IS_PRESENT_IN_DATASHEET}{...}
#'   \item{IS_PRESENT_BUT_NOT_IN_DS}{...}
#'   \item{VIAL_BARCODE}{...}
#'   \item{CONTAINER_ARRAY_TYPE_ID}{...}
#'   \item{DNA_TRAY_WORKBENCH_ID}{...}
#'   \item{DNA_TRAY_CODE}{...}
#'   \item{DNA_TRAY_WELL_POS}{...}
#'   \item{DNA_TRAY_WELL_CODE}{...}
#'   \item{STORAGE_ID}{...}
#'   \item{UNIT}{...}
#'   \item{SHELF_RACK}{...}
#'   \item{SLOT}{...}
#'   \item{EXHAUSTED_HOW}{...}
#'   \item{EXHAUSTED_BY}{...}
#'   \item{EXHAUSTED_DATE}{...}
#'   \item{AGENCY}{...}
#'   \item{OTHER_AGENCY_KEY}{...}
#'   \item{NUM_OTOLITHS_MISSING}{...}
#'   \item{OTO_INVENTORY_COMMENT}{...}
#' }
#'
#' If \code{import.vars = FALSE}, the function returns a tibble with 7 collection information variables and the 31 tissue import variables.
#' The additional collection information variables are:
#' \itemize{
#'   \item{SILLY_CODE}{...}
#'   \item{REGION_CODE}{...}
#'   \item{QUADRANT}{...}
#'   \item{LOCATION_CODE}{...}
#'   \item{LOCATION_DESCRIPTOR}{...}
#'   \item{LIFE_STAGE}{...}
#'   \item{COLLECTION_TYPE}{...}
#' }
#'
#' @note
#' Output of this function can be used for creating an import file for the Loki tissue importer, 
#' adding additional information to the Loki tissue table (default output), and 
#' generating extraction lists for lab projects.
#'
#' The tissue importer requires that all 31 columns are in the correct order and are present. 
#' When writing the import file, ensure it does not contain NAs. To export the tissue import file 
#' using [readr::write_csv()], set the `eol` argument to "\\r\\n" for proper formatting:
#' \code{import_file %>% readr::write_csv(file = "ImportFile.csv", na = "", eol = "\\r\\n")}
#'
#' Additionally, if you have date columns, change their format to 'numeric' in Microsoft Excel 
#' before importing tissue information with dates.
#' 
#' @details
#' There are several ways to pull tissue data with this function. You can pull by sillyvec, by unit, or by unit and shelf.rack. 
#' Any other combination of those arguments will throw an error message.
#' 
#'
#' @examples
#' \dontrun{
#' #By sillyvec
#' get_tissue_data(
#'   sillyvec = c("KCDVF18", "KCDVF19", "KCDVF20"),
#'   unit = NULL,
#'   shelf.rack = NULL,
#'   username = username,
#'   password = password,
#'   file = path.expand("~/test_tissue_table.csv"),
#'   import.vars = FALSE
#' )
#' 
#' #By unit
#'  get_tissue_data(
#'   sillyvec = NULL,
#'   unit = "WQ",
#'   shelf.rack = NULL,
#'   username = username,
#'   password = password,
#'   file = path.expand("~/test_tissue_table.csv"),
#'   import.vars = FALSE
#' )
#' 
#' #By unit and shelf.rack
#'  get_tissue_data(
#'   sillyvec = NULL,
#'   unit = "WQ",
#'   shelf.rack = 1:10,
#'   username = username,
#'   password = password,
#'   file = path.expand("~/test_tissue_table.csv"),
#'   import.vars = FALSE
#' )
#' 
#' }
#'
#' @export
get_tissue_data <- function(sillyvec = NULL, unit = NULL, shelf.rack = NULL, username, password, file = NULL, import.vars = TRUE) {
  
  if(!is.null(sillyvec) & !is.null(unit) & !is.null(shelf.rack)| 
     !is.null(sillyvec) & is.null(unit) & !is.null(shelf.rack)|
     !is.null(sillyvec) & !is.null(unit) & is.null(shelf.rack)|
     is.null(sillyvec) & is.null(unit) & !is.null(shelf.rack)){
    
    stop("Only supply one of the following: sillyvec, unit, or unit and shelf.rack")
    
  }
  
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  #By sillyvec
  if(!is.null(sillyvec) & is.null(unit) & is.null(shelf.rack)){
    
    qry <- paste("SELECT * FROM AKFINADM.V_GEN_SAMPLED_FISH_TISSUE WHERE SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")", sep = "") #Query
    
  }
  
  #By unit

  if(is.null(sillyvec) & !is.null(unit) & is.null(shelf.rack)){
    
    qry <- paste("SELECT * FROM AKFINADM.V_GEN_SAMPLED_FISH_TISSUE WHERE UNIT IN (", paste0("'", unit, "'", collapse = ","), ")", sep = "") #Query
    
  }
  
  #By unit and self.rack
  
  if(is.null(sillyvec) & !is.null(unit) & !is.null(shelf.rack)){
    
    qry <- paste("SELECT * FROM AKFINADM.V_GEN_SAMPLED_FISH_TISSUE WHERE UNIT IN (", paste0("'", unit, "'", collapse = ",") , ") AND SHELF_RACK IN (", paste0("'", shelf.rack, "'", collapse = ","), ")", sep = "") #Query
    
  }
  
  dataAll0 <- RJDBC::dbGetQuery(con, qry)  #Pulling data from LOKI using the connection and tissue data query
  
  discon <- RJDBC::dbDisconnect(con) # Disconnect from Loki
  
  #Replace 0 with NA for variables that have a check box in Loki. If zeros are present in the import they get treated the same as 1's
  dataAll <- dataAll0 %>% 
    dplyr::mutate(dplyr::across(dplyr::all_of(c("IS_MISSING_PAIRED_DATA_EXISTS", "IS_PRESENT_IN_DATASHEET", "WELL_HAS_MORE_THAN_ONE_SAMPLE", "IS_PRESENT_BUT_NOT_IN_DS")) , ~gsub(pattern = 0, replacement = NA, x = .)))
  
  if (import.vars == TRUE) {
    
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
  
  if (!is.null(file)) {
    
    readr::write_csv(output, file, na = "") # Write out a csv file without NAs
    
  }
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime)
  
  return(output)
}