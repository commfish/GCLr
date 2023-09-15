#' Get a GTSeq Project Sample Sheet
#' 
#' This function connects to Loki and gets the GTSeq sample sheet for a given project. 
#' i7 barcodes are added to the sheet starting with the supplied `i7_NAME_start`. (see details)
#' 
#' @param project the project name from Loki iStrategy (e.g., "P015")
#' @param i7_NAME_start a character string with three characters (i.e., `nchar(i7_NAME_start) = 3`) and leading zeros.
#' @param i7_NAME_max this is the maximum i7 name to use when creating the sample sheet. 
#'                    This may change periodically depending on what barcodes the lab has in stock. Allowable values are `1:150`.  
#'                    Check with lab staff if you aren't sure what to put here. 
#' @param dir the directory (folder) where you want to write out the sample sheet .csv file. \emph{Note: The file will be named automatically using the supplied project name.}
#' @param username your state user name
#' @param password your password used to access LOKI; see Eric Lardizabal if you need to set up a password
#' 
#' @details
#' Before running this function, the project information must be entered into `iStrategy`. 
#' 
#' i7 barcodes come in individual vials and get used up over time. 
#' To make sure the vials are used evenly (i.e., run out at the same time), each project will start with the next i7_NAME after the last i7 name that was use in the previous project. 
#' To determine which what i7 barcode that was used last, look in "S:/GCL/Lab/Genotyping/SNP Projects/NextSeq Flow Cell Log.xlsx".
#' The last project should be at the bottom of the sheet if the "Run Date" field is sorted lowest to highest. The "i7s" field will contain the range of i7s used in the project.   
#' 
#' @note as of 9/12/23 the lab is only using i7 barcodes '001'-'075'.
#' 
#' @returns This function produces a tibble containing the following fields:
#'    \itemize{
#'                \item \strong{Project}: project name
#'                \item \strong{SILLY_FISH_ID}: silly code and fish number
#'                \item \strong{Plate_ID}: extraction plate ID number
#'                \item \strong{i7name}: i7 barcode name for each plate
#'                \item \strong{i7sequence}: i7 barcode sequence
#'                \item \strong{WellID}: extraction plate well ID (e.g., A1, B1, C1, ....., H12)
#'                \item \strong{i5name}: i5 barcode name for each fish
#'                \item \strong{i5sequence}: i5 barcode sequence
#'                \item \strong{Panel}: the marker panel name (marker suite)
#'                }
#'                
#'          The tibble is then written out to a .csv file.
#' 
#' @examples
#' \dontrun{
#' 
#'  get_gtseq_sample_sheet (project = "P015", i7_NAME_start = "001", i7_NAME_max = 75, dir = path.expand("~"), username = "awbarclay", password = scan(path.expand("~/R/usr.pw"), what = "")[2])
#' 
#' }
#'            
#' @export
get_gtseq_sample_sheet <- function(project, i7_NAME_start, i7_NAME_max = 75, dir, username, password){
  
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  procedure <- paste0("BEGIN AKFINADM.POPULATE_TEMP_GTSEQ_SAMPLE_SHEET_CSV('", project, "', '", i7_NAME_start, "', ", as.numeric(i7_NAME_max), ", '", username, "'); END;") #Sample sheet procedure
  
  #procedure <- paste0("BEGIN AKFINADM.POPULATE_TEMP_GTSEQ_SAMPLE_SHEET_CSV('", project, "', '", i7_NAME_start, "', ", as.numeric(i7_NAME_max), ", 'blah'); END;") #Sample sheet procedure
  
  RJDBC::dbSendUpdate(con, procedure)  #Pulling data from LOKI using the connection and genotype query
  
  my.table <- RJDBC::dbGetQuery(con, paste0("SELECT * FROM AKFINADM.TEMP_GTSEQ_SAMPLE_SHEET_CSV_", username, " ORDER BY PLATE_ID, I5_NAME"))
  
  #my.table <- RJDBC::dbGetQuery(con, paste0("SELECT * FROM AKFINADM.TEMP_GTSEQ_SAMPLE_SHEET_CSV_blah ORDER BY PLATE_ID, I5_NAME"))
  
  discon <- RJDBC::dbDisconnect(con)
  
  if(nrow(my.table)==0){

    stop("No inforation exists for project ", project, ". Check to make sure the spelling is correct and that the project iformation has been entered in iStrategy.")

  }

  file <- paste0(dir, "/GTSeq_SampleSheet_FlowcellID_", project, ".csv")
  
  readr::write_csv(my.table, file = file)
  
  return(my.table)
  
  print(Sys.time()-start.time)
  
}  