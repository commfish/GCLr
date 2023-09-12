#' Get a GTSeq Project Sample Sheet
#' 
#' This function connects to Loki and gets the GTSeq sample sheet for a given project. 
#' i7 barcodes are added to the sheet starting with the supplied `i7_NAME_start`. (see details)
#' 
#' @param project the project name from Loki iStrategy (e.g., "P015")
#' @param i7_NAME_start a character string with three characters (i.e., `nchar(i7_NAME_start) = 3`) and leading zeros.  
#' @param dir the directory (folder) where you want to write out the sample sheet .csv file. \emph{Note: The file will be named automatically using the supplied project name.}
#' @param username your state user name
#' @param password your password used to access LOKI; see Eric Lardizabal if you need to set up a password
#' 
#' @details
#' For this function to work, the project information must be entered into [iStrategy] first. 
#' i7 barcodes come in individual vials and get used up over time. 
#' To make sure the vials are used evenly (i.e., run out at the same time), each project will start with the next i7_NAME after the last i7 name that was use in the previous project. 
#' To determine which what i7 barcode that was used last, look 'NextSeq Flow Cell Log.xlsx' located here ["S:\GCL\Lab\Genotyping\SNP Projects"]("S:/GCL/Lab/Genotyping/SNP Projects").
#' The last project should be at the bottom of the sheet if the "Run Date" field is sorted lowest to highest. The "i7s" field will contain the range of i7s used in the project.   
#' 
#' @returns This function writes out a .csv file containing the following fields:
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
#' @examples
#' \dontrun{
#' 
#'   get_gtseq_sample_sheet (project = "P015", i7_NAME_start = "001", dir = path.expand("~"), username = "awbarclay", password = scan(path.expand("~/R/usr.pw"), what = "")[2])
#' 
#' }
#'            
#' @export
get_gtseq_sample_sheet <- function(project, i7_NAME_start, dir, username, password){
  
  if(!is.character(i7_NAME_start)| !nchar(i7_NAME_start)==3){
    
    stop("The 'i7_NAME_start' supplied needs to be a character string with three characters and leading zeros. (e.g. '001')")
    
  }
  
  #This message will need to be updated if the lab ever gets barcodes beyond '075'.
  if(as.numeric(i7_NAME_start) > 75){
    
    stop("Currently, the only allowable i7 barcodes are '001'-'075'. Please supply a i7 name within that range and try again. 
         If the lab has i7 barcodes beyond '075', you will need to send a message to the labs analyst programer to add the additional barcodes to the 
         table in Loki and start an issue on the GCLr repository 'https://github.com/commfish/GCLr' to update this function.") 
    
  }
  
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  procedure <- paste0("BEGIN AKFINADM.POPULATE_TEMP_GTSEQ_SAMPLE_SHEET_CSV('", project, "', '", i7_NAME_start, "'); END;") #Sample sheet procedure
  
  RJDBC::dbSendUpdate(con, procedure)  #Pulling data from LOKI using the connection and genotype query
  
  my.table <- RJDBC::dbGetQuery(con, "SELECT * FROM AKFINADM.TEMP_GTSEQ_SAMPLE_SHEET_CSV ORDER BY PLATE_ID, I5_NAME")
  
  discon <- RJDBC::dbDisconnect(con)
  
  if(nrow(my.table)==0){
    
    stop("No inforation exists for project ", project, ". Check to make sure the spelling is correct and that the project iformation has been entered in iStrategy.")
    
  }
  
  file <- paste0(dir, "/GTSeq_SampleSheet_FlowcellID_", project, ".csv")
  
  readr::write_csv(my.table, file = file)
  
  print(Sys.time()-start.time)
  
}  