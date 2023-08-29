#' @title Create GCL Objects for a Lab Project
#'
#' @description This function is intended for use in the qc.R script. It pulls project genotypes to create "slim" .gcl objects for each silly and the LocusControl for all loci used in the project.
#'              
#' @param project_name A character vector of one or more project names as spelled in LOKI.
#' @param sillyvec A character vector of SILLYs.
#' @param loci A character vector of locus names as they are spelled in LOKI.
#' @param username Your username for accessing LOKI through R.
#' @param password Your password for accessing LOKI through R.
#'
#' @return Returns `project_sillys`, a character vector of all sillys in the project (sillyvec).
#' @return Returns "slim" .gcl objects for each silly (slim = not all attributes table).
#' @return Returns `LocusControl` for all loci used in the project (including "loci", "nallales", "ploidy", "alleles").
#'
#' @details
#' This function pulls project genotypes from LOKI database to create "slim" .gcl objects for each silly, along with the LocusControl for all loci used in the project.
#' The function checks the combination of arguments to ensure the correct usage.
#' 
#' @note
#' Warning: Genotypes (pre-October 2016) may not exist in the LOKI lookup table with a project name, so genotypes will have to be pulled by sillyvec and loci.
#' 
#' @examples
#' \dontrun{
#' read_proj_geno(project_name = c("P014", "P015", "P016"), sillyvec = c("SCIMA18", "SCIMA17"), 
#'                loci = c("One_E2", "One_MHC2_251", "One_Cytb_17"), username = "awbarclay", password = "password")
#' }
#'
#' @export
loki2r_proj <- function(project_name = NULL, sillyvec = NULL, loci = NULL, username, password) {
 
    # Recording function start time
  start.time <- Sys.time()
  
  if (exists("LocusControl", where = 1)) {
    
    stop("LocusControl already exists")
    
  }
  
  # Checking to make sure the correct combination of arguments is being use. 
  # If wrong, the function will stop and print an error message to the console.
  if (is.null(sillyvec) & is.null(loci) & is.null(project_name) |
      !is.null(sillyvec) & !is.null(loci) & !is.null(project_name) |
      !is.null(sillyvec) & is.null(loci) & !is.null(project_name) |
      is.null(sillyvec) & !is.null(loci) & !is.null(project_name) |
      is.null(sillyvec) & !is.null(loci) & is.null(project_name))
  {
    stop("The user must supply one of the following argument combinations:\n  1) sillyvec (for all loci and individuals for each silly),\n  2) sillyvec and loci (all individuals for supplied locus list), or\n  3) project_name (for all individuals and loci in a given project)")
  }
  
  # Setting default java.parameters
  options(java.parameters = "-Xmx10g")
  
  # Connect to LOKI
  url <- GCLr:::loki_url() #Function that gets the correct URL to access the database on the oracle cloud
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  

  # Get genotypes
  # Creating java query when sillyvec and loci are supplied.  
  if (!is.null(sillyvec) & !is.null(loci)) {
    gnoqry <- paste0("SELECT * FROM AKFINADM.R_READ_PROJECT_GENOTYPES WHERE LOCUS IN (", paste0("'", loci, "'", collapse = ","), ") AND SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")")
  } 
  
  # Creating java query when only sillyvec is supplied
  if (!is.null(sillyvec) & is.null(loci)) {
    gnoqry <- paste0("SELECT * FROM AKFINADM.R_READ_PROJECT_GENOTYPES WHERE SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")")
  }
  
  # Creating java query when only project_name is supplied.
  if (!is.null(project_name)) {
    gnoqry <- paste0("SELECT * FROM AKFINADM.R_READ_PROJECT_GENOTYPES GENO WHERE EXISTS (SELECT * FROM AKFINADM.V_LAB_PROJECT_WELL LPW WHERE LPW.LAB_PROJECT_NAME IN (", paste0("'", project_name, "'", collapse = ","), ") AND LPW.SILLY_CODE = GENO.SILLY_CODE AND LPW.FISH_NO = GENO.FK_FISH_ID)")
  }
  
  # Pull genotypes
  dataAll <- RJDBC::dbGetQuery(con, gnoqry) %>% 
    dplyr::as_tibble()
  
  # Filter for project_name again, if applicable
  if (!is.null(project_name)) {
    dataAll <- dataAll %>% 
      dplyr::filter(LAB_PROJECT_NAME == project_name)
  }
  
  # Get list of unique sillys and assign `project_sillys`; this is needed for qc script
  sillyvec <- unique(dataAll$SILLY_CODE)
  assign(x = "project_sillys", value = sillyvec, pos = 1)
  
  # Disconnect from LOKI
  discon <- RJDBC::dbDisconnect(con)
  

  # Get list of unique loci and assign `LocusControl`; this is needed for qc script
  loci <- sort(unique(dataAll$LOCUS))
  
  GCLr::create_locuscontrol(locusnames = loci, username = username, password = password)
  
# This section generates the .gcl objects  
  # xxx CSJ - shouldn't I be able to pivot without the left join??
  # xxx CSJ - make sure this is correct,paying extra attention to handling 0s or NAs in Alleles
  # xxx CSJ - why does that note below (original code) say this still needs work? Maybe it got left by accident, or perhaps if you don't use the project_name pull?
  
  #convert data to wide format tibble
  data_master <- dataAll %>% 
    tidyr::unite(SillySource, SILLY_CODE, FK_FISH_ID, sep = "_", remove = FALSE) %>% 
    dplyr::select(-ALLELE_1, -ALLELE_2) %>%
    dplyr::rename(ALLELE_1 = ALLELE_1_FIXED, ALLELE_2 = ALLELE_2_FIXED) %>% 
    dplyr::mutate(ALLELE_1 = gsub(pattern = "0", replacement = NA_character_, x = ALLELE_1),
                  ALLELE_2 = gsub(pattern = "0", replacement = NA_character_, x = ALLELE_2)) %>% 
    tidyr::pivot_wider(id_cols = SillySource, names_from = LOCUS, values_from = c(ALLELE_1, ALLELE_2), names_glue = "{LOCUS}_{.value}", values_fill = NULL)  %>% #note values_fill to keep na as na
    dplyr::rename_with(~ gsub("_ALLELE_1", "", .), ends_with("_ALLELE_1")) %>% #fix names locus
    dplyr::rename_with(~ gsub("_ALLELE_2", ".1", .), ends_with("_ALLELE_2")) %>% #fix names locus.1
    dplyr::left_join(dplyr::distinct(dataAll %>% 
                       dplyr::select(c("FK_FISH_ID", "COLLECTION_ID", "SILLY_CODE", "PLATE_ID", "PK_TISSUE_TYPE")) %>% 
                       tidyr::unite(SillySource, SILLY_CODE, FK_FISH_ID, sep = "_", remove = FALSE)),
                     by = dplyr::join_by(SillySource)) #join back in missing columns from pivot
    
  # Add extra columns as NAs for compatibility with all other functions
  attnames <- c("FK_FISH_ID", "COLLECTION_ID", "SILLY_CODE", "PLATE_ID", "PK_TISSUE_TYPE", 
                "CAPTURE_LOCATION", "CAPTURE_DATE", "END_CAPTURE_DATE", "MESH_SIZE", 
                "MESH_SIZE_COMMENT", "LATITUDE", "LONGITUDE", "AGENCY", "VIAL_BARCODE", 
                "DNA_TRAY_CODE", "DNA_TRAY_WELL_CODE", "DNA_TRAY_WELL_POS", "CONTAINER_ARRAY_TYPE_ID", 
                "SillySource") # list of 19 allowed attributes
  
  # Check which attnames exist in the data_master tibble, identify missing
  missing_columns <- setdiff(attnames, names(data_master))
  
  # Add the missing attnames, as NA values, and arrange (19 attributes first)
  data_master <- data_master %>%
    tibble::add_column(!!!setNames(rep(NA, length(missing_columns)), missing_columns)) %>% 
    dplyr::select(tidyselect::all_of(attnames), sort(names(.), na.last = TRUE))
    
    
  message("Data successfully pulled from LOKI, building SILLY.gcl objects.\n")
  
  message("Note - this is a slim gcl object, with all 19 attributes added (as NA) for compatibility with other functions.\n")
  
  # Make .gcl objects by silly
  for (silly in sillyvec) {  
    
    data_silly <- dplyr::filter(data_master, SILLY_CODE == silly)
    
    # Adaptation of LOKI2R to allow pulling all SILLYs just based on ProjectID
    if (nrow(data_silly) == 0) {
      warning(paste(silly, "was not found in LOKI, it may not have been imported!!!", sep = " "))
      message(paste0(silly,".gcl NOT FOUND IN LOKI ", match(silly,sillyvec)," of ",length(sillyvec)," FAILED!!!"))
      next
    }
  
    ids <- as.character(sort(unique(data_silly$FK_FISH_ID)))
    nind <- length(ids)
    
    if (!nind) {
      message0 <- paste0(silly," is empty.")
      message(message0)
      next
    }
    
    message(paste0(silly, ".gcl created ", match(silly, sillyvec), " of ", length(sillyvec), " completed."))
    
    assign(x = paste0(silly, ".gcl"),value = data_silly, pos = 1)
   
  }  # silly
  
  stop.time <- Sys.time()
  fulltime <- stop.time - start.time
  print(fulltime) 
}