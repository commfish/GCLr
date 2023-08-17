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
#' @return Returns `ProjectSillys`, a character vector of all sillys in the project (sillyvec).
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
  
  # Get list of unique sillys and assign `ProjectSillys`; this is needed for qc script
  sillyvec <- unique(dataAll$SILLY_CODE)
  assign(x = "ProjectSillys", value = sillyvec, pos = 1)
  
  # Disconnect from LOKI
  discon <- RJDBC::dbDisconnect(con)
  

  # Get list of unique loci and assign `LocusControl`; this is needed for qc script
  loci <- sort(unique(dataAll$LOCUS))
  #nloci <- length(loci) # xxx CSJ I don't think this is needed anymore
  
  GCLr::create_locuscontrol(locusnames = loci, username = username, password = password)
  
  assign(x = "loci", value = LocusControl$locusnames, pos = 1)
  assign(x = "nalleles", value = LocusControl$nalleles, pos = 1)
  assign(x = "ploidy", value = LocusControl$ploidy, pos = 1)
  assign(x = "alleles", value = LocusControl$alleles, pos = 1)
  
  
# This section generates the .gcl objects  
  # xxx CSJ - shouldn't I be able to pivot without the left join??
  # xxx CSJ - make sure this is correct,paying extra attention to handling 0s or NAs in Alleles
  # xxx CSJ - why does that note below (original code) say this still needs work? Maybe it got left by accident, or perhaps if you don't use the project_name pull?
  
  #convert data to wide format tibble
  data_master <- dataAll %>% 
    tidyr::unite(SillySource, SILLY_CODE, FK_FISH_ID, sep = "_", remove = FALSE) %>% 
    dplyr::select(-ALLELE_1, -ALLELE_2) %>%
    dplyr::rename(ALLELE_1 = ALLELE_1_FIXED, ALLELE_2 = ALLELE_2_FIXED) %>% 
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
  
###===============   
  
#   
#   
#   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # Join genotypes with extraction information (PLATE_ID) and make .gcl objects
#   # This still needs work because it pulls all plates per fish, not just the project plates!!!
#   
#   attnames <- c("FK_FISH_ID", "COLLECTION_ID", "SILLY_CODE", "PLATE_ID", "PK_TISSUE_TYPE", "SillySource")
#   
#   data_master <- dataAll %>% 
#     tidyr::unite(SillySource, SILLY_CODE, FK_FISH_ID, sep = "_", remove = FALSE) %>% 
#     dplyr::select(-ALLELE_1, -ALLELE_2) %>%
#     dplyr::rename(ALLELE_1 = ALLELE_1_FIXED, ALLELE_2 = ALLELE_2_FIXED) %>% 
#     tidyr::unite(GENO, ALLELE_1, ALLELE_2, sep = "/", remove = FALSE) %>% 
#     dplyr::mutate(ALLELES = dplyr::case_when(PLOIDY == "D" ~ GENO,
#                                              PLOIDY == "H" ~ ALLELE_1)) %>% 
#     dplyr::select(c(dplyr::all_of(attnames), "LOCUS", "ALLELE_1", "ALLELE_2", "ALLELES")) 
#   
#   message("Data successfully pulled from LOKI, building SILLY.gcl objects\n")
#   
#   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # Make .gcl objects by silly
#   for (silly in sillyvec) {  
#     
#     data_silly <- dplyr::filter(data_master, SILLY_CODE == silly)
#     
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # Adaptation of LOKI2R to allow pulling all SILLYs just based on ProjectID
#     if (nrow(data_silly) == 0) {
#       warning(paste(silly, "was not found in LOKI, it may not have been imported!!!", sep = " "))
#       message(paste0(silly,".gcl NOT FOUND IN LOKI ", match(silly,sillyvec)," of ",length(sillyvec)," FAILED!!!"))
#       next
#     }
#     
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     
#     ids <- as.character(sort(unique(data_silly$FK_FISH_ID)))
#     nind <- length(ids)
#     
#     if (!nind) {
#       message0 <- paste0(silly," is empty.")
#       message(message0)
#       next
#     }
#     
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # Build scores array
#     
#     scores <- array(data = NA,
#                     dim = c(nind, nloci, max(ploidy)), 
#                     dimnames = list(ids, loci, paste0("Dose", seq(max(ploidy))))
#     )
#     
#     sillyloci <- sort(unique(data_silly$LOCUS))  # These are the loci available for the current silly, this is needed to subset the tapply
#     
#     scores_allele_1 <- data_silly %>% 
#       dplyr::filter(LOCUS %in% sillyloci) %>% 
#       dplyr::select(FK_FISH_ID, LOCUS, ALLELE_1) %>% 
#       dplyr::mutate(FK_FISH_ID = as.character(FK_FISH_ID)) %>% 
#       tidyr::spread(LOCUS, ALLELE_1) %>% 
#       tibble::column_to_rownames(var = "FK_FISH_ID") %>% 
#       as.matrix()
#     
#     scores[ids, sillyloci, "Dose1"] <- scores_allele_1[ids, sillyloci]
#     
#     scores_allele_2 <- data_silly %>% 
#       dplyr::filter(LOCUS %in% sillyloci) %>% 
#       dplyr::select(FK_FISH_ID, LOCUS, ALLELE_2) %>% 
#       dplyr::mutate(FK_FISH_ID = as.character(FK_FISH_ID)) %>% 
#       tidyr::spread(LOCUS, ALLELE_2) %>% 
#       tibble::column_to_rownames(var = "FK_FISH_ID") %>% 
#       as.matrix()
#     
#     scores[ids, sillyloci, "Dose2"] <- scores_allele_2[ids, sillyloci]
#     
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # Build counts array
#     
#     counts <- array(data = NA, 
#                     dim = c(nind, nloci, max(nalleles)), 
#                     dimnames = list(ids, loci, paste0("Allele ", seq(max(nalleles))))
#     )
#     
#     for (ind in ids) {
#       for (locus in loci) {
#         for (al in seq(nalleles[locus])) {
#           counts[ind, locus, al] <- sum(scores[ind, locus, seq(ploidy[locus])] == alleles[[locus]][al])
#         }  # al
#       }  # locus
#       counts0 = counts[ind, , ]
#       counts0[is.na(counts0)] <- 0
#       zeroBOOL <- apply(counts0,1,sum) != ploidy
#       counts[ind, zeroBOOL, ] <- NA
#     }  # ind
#     
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # Build attributes data.frame
#     
#     attributes <- data_silly %>% 
#       dplyr::select(dplyr::all_of(attnames)) %>% 
#       dplyr::distinct() %>% 
#       dplyr::mutate(FISH_ID = as.character(FK_FISH_ID),
#                     PLATE_ID = as.character(PLATE_ID)) %>% 
#       dplyr::arrange(FK_FISH_ID) %>% 
#       tibble::column_to_rownames(var = "FISH_ID")
#     
#     message(paste0(silly, ".gcl created ", match(silly, sillyvec), " of ", length(sillyvec), " completed."))
#     
#     assign(paste0(silly, ".gcl"), list(counts = counts, scores = scores, n = nind, attributes = attributes), pos = 1)
#     
#   }  # silly
#   
#   stop.time <- Sys.time()
#   fulltime <- stop.time - start.time
#   print(fulltime) 
# }