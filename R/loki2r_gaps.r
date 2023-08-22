#' Create GCL Objects (GAPS)
#'
#' This function connects to LOKI and creates a ".gcl" object for each silly code in `sillyvec`. A "*.gcl" object is a Gene Conservation Laboratory genotypes object with associated sample attributes.
#'
#' @param sillyvec Vector of silly codes to pull from LOKI.
#' @param username Username for LOKI connection.
#' @param password Password for LOKI connection.
#' 
#' @return The runtime of the function in seconds 
#' 
#' Assigns the following objects to the global environment:
#'  
#' A tibble of individuals missing data at one or more loci called `missing_indvs_loci`. 
#' 
#' A tibble with the following columns for each silly:
#'    \itemize{
#'                \item \code{FK_FISH_ID} (double): fish ID numbers for each individual
#'                \item \code{COLLECTION_ID} (double): the unique collection ID number for each individual
#'                \item \code{SILLY_CODE} (character): the silly code for each individual
#'                \item \code{PLATE_ID} (character): the extraction plate ID for each individual
#'                \item \code{PK_TISSUE_TYPE} (character): the tissue type extracted for DNA for each individual
#'                \item \code{CAPTURE_LOCATION} (character): the location where each individual was captured for sampling
#'                \item \code{CAPTURE_DATE} (date): the date each individual was captured (e.g., May 5, 2020 = "2020-05-05")
#'                \item \code{END_CAPTURE_DATE} (date): the last collection date for a silly (e.g., May 5, 2020 = "2020-05-05")
#'                \item \code{MESH_SIZE} (character): the mesh size of the net used to capture (harvest) each individual
#'                \item \code{MESH_SIZE_COMMENT} (character): comments about mesh size
#'                \item \code{LATITUDE} (double): the latitude where each individual was captured in decimal degrees
#'                \item \code{LONGITUDE} (double): the longitude where each individual was captured in decimal degrees
#'                \item \code{AGENCY} (character): the name of the agency or organization that collected each individual
#'                \item \code{VIAL_BARCODE} (character): the barcode on the collection vial
#'                \item \code{DNA_TRAY_CODE} (character): the barcode on the collection tray/card
#'                \item \code{DNA_TRAY_WELL_CODE} (double): the unique number assigned to each position in the collection tray/card for each individual (e.g., positions A1-A10 = codes 1-10)
#'                \item \code{DNA_TRAY_WELL_POS} (character): the position in the collection tray/card (e.g., A1, A2, B1, B2, etc.)
#'                \item \code{CONTAINER_ARRAY_TYPE_ID} (double): the number code for the collection container (e.g., tray or card)
#'                \item \code{SillySource} (double): the original silly code and fish ID for each individual (e.g., KQUART06_1). When pulled from loki this will be the SILLY_CODE and FK_FISH_ID
#'                \item \code{genotypes} with a column for each allele for each locus
#'                \item \code{CAPTURE_LOCATION} (character): the location where each individual was captured for sampling
#'                \item \code{CAPTURE_DATE} (date): the date each individual was captured (e.g., May 5, 2020 = "2020-05-05")
#'                \item \code{END_CAPTURE_DATE} (date): the last collection date for a silly (e.g., May 5, 2020 = "2020-05-05")
#'                \item \code{MESH_SIZE} (character): the mesh size of the net used to capture (harvest) each individual
#'                \item \code{MESH_SIZE_COMMENT} (character): comments about mesh size
#'                \item \code{LATITUDE} (double): the latitude where each individual was captured in decimal degrees
#'                \item \code{LONGITUDE} (double): the longitude where each individual was captured in decimal degrees
#'                \item \code{AGENCY} (character): the name of the agency or organization that collected each individual
#'                \item \code{VIAL_BARCODE} (character): the barcode on the collection vial
#'                \item \code{DNA_TRAY_CODE} (character): the barcode on the collection tray/card
#'                \item \code{DNA_TRAY_WELL_CODE} (double): the unique number assigned to each position in the collection tray/card for each individual (e.g., positions A1-A10 = codes 1-10)
#'                \item \code{DNA_TRAY_WELL_POS} (character): the position in the collection tray/card (e.g., A1, A2, B1, B2, etc.)
#'                \item \code{CONTAINER_ARRAY_TYPE_ID} (double): the number code for the collection container (e.g., tray or card)
#'                \item \code{SillySource} (double): the original silly code and fish ID for each individual (e.g., KQUART06_1). When pulled from loki this will be the SILLY_CODE and FK_FISH_ID
#'                \item \code{genotypes} with a column for each allele for each locus
#'              }
#'                
#' @details This function requires the R package "RJDBC" for connecting to LOKI.
#' It connects to LOKI using the provided \code{username} and \code{password}, and retrieves genotypes and sample attributes for each silly code in \code{sillyvec}.
#' The genotypes and attributes are stored in separate "*.gcl" objects named after each silly code.
#'
#' @examples
#' \dontrun{
#' sillyvec <- c("KQUART06", "KCRESC06", "KJUNE05")
#' loki2r_gaps(sillyvec = sillyvec, username = "awbarclay", password)
#' }
#' 
#' @export
loki2r_gaps <- function(sillyvec, username, password){

  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  markersuite <- "GAPS_Chinook_uSATs"
  
  lociqry <- paste0("SELECT * FROM AKFINADM.V_LOCUS_MARKERSUITE WHERE SUITE_NAME = '", markersuite, "'")
  
  locidata <- RJDBC::dbGetQuery(con, lociqry)
  
  locusnames <- as.character(locidata$LOCUS_NAME)
  
  Publishedlocusnames <- as.character(locidata$PUBLISHED_NAME)
  
  nloci <- length(locusnames)
  
  ploidy <- sapply(locusnames, function(locus) {
    qry <- paste0("SELECT PLOIDY, LOCUS_NAME FROM AKFINADM.LOCUS_LOOKUP WHERE LOCUS_NAME IN ", "'", locus, "'")
    2^(as.character(RJDBC::dbGetQuery(con, qry)$PLOIDY[1]) == "D")
  } )
  
  ConTable <- RJDBC::dbGetQuery(con, "SELECT * FROM AKFINADM.GAPS_ALLELE_CONVERSION")
  
  ConTable$VALUE_ADFG <- as.numeric(ConTable$VALUE_ADFG)  # this is crucial for allele sorting so that "99" is at the front, not the back of the allele list
  
  ConTable$VALUE_CTC <- as.numeric(ConTable$VALUE_CTC)  # this is crucial for allele sorting so that "99" is at the front, not the back of the allele list
  
  alleles <- sapply(locusnames, function(loc) {
    
    as.character(sort(unique(as.vector(ConTable[ConTable[, "LOCUS_NAME"] == loc, "VALUE_CTC"]))))
    
    })
  
  nalleles <- sapply(alleles, function(allele) {length(allele)} )
  
  LocusCtl <- tibble::tibble(MarkerSuite = markersuite, locusnames = locusnames, Publishedlocusnames = Publishedlocusnames, nalleles = nalleles, ploidy = ploidy, alleles = alleles)
  
  assign("LocusControl", LocusCtl, pos = 1, envir = .GlobalEnv) #Assign elements to LocusControl list.	
  
  loci <- LocusCtl$locusnames
  
  nloci <- length(loci)
  
  gnoqry.norm <- paste("SELECT * FROM AKFINADM.V_GNOQRY WHERE LOCUS IN (", paste0("'", loci, "'", collapse = ","), ") AND SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")", sep = "") #Gentoype query
  
  dataAll.norm <- RJDBC::dbGetQuery(con, gnoqry.norm) %>% 
    tibble::as_tibble() 
  
  gnoqry.gaps <- paste0("SELECT * FROM AKFINADM.V_GAPS_GENOTYPES WHERE SUITE_NAME = '", markersuite, "' AND FK_COLLECTION_ID IN (", paste0("'", unique(dataAll.norm$COLLECTION_ID), "'", collapse = ","), ") ORDER BY LOCUS, FK_COLLECTION_ID,FK_FISH_ID")
  
  dataAll.gaps <- RJDBC::dbGetQuery(con, gnoqry.gaps) %>% 
    tibble::as_tibble() %>% 
    dplyr::select(LOCUS, FK_COLLECTION_ID, FK_FISH_ID, ALLELE1_CONV, ALLELE2_CONV)
  
  dataAll0 <- dplyr::inner_join(dataAll.norm, dataAll.gaps, by = c("COLLECTION_ID" = "FK_COLLECTION_ID", "LOCUS", "FISH_ID"="FK_FISH_ID")) %>% 
    dplyr::mutate(ALLELE_1 = ALLELE1_CONV, ALLELE_2 = ALLELE2_CONV) %>% 
    dplyr::select(-ALLELE1_CONV, -ALLELE2_CONV)
  
  dataAllbool <- dataAll0$PK_TISSUE_TYPE == dataAll0$PREFERRED_TISSUE 
  
  bothNAbool <- is.na(dataAll0$PREFERRED_TISSUE) &  is.na(dataAll0$PK_TISSUE_TYPE)
  
  dataAllbool[is.na(dataAllbool)] <- FALSE 
  
  dataAllbool[bothNAbool] <- TRUE  
  
  dataAll <- dataAll0[dataAllbool, ] %>% 
    tibble::as_tibble() %>% 
    dplyr::select(-PREFERRED_TISSUE)
  
  discon <- RJDBC::dbDisconnect(con)
  
  missing_sillys <- setdiff(sillyvec, dataAll$SILLY_CODE %>% unique()) #Find which sillys had no data for any loci in LocusControl
  
  # what indvs are missing loci from LocusControl (i.e. no genotyping attempted)?
  missing_indvs <- dataAll %>% 
    dplyr::count(SILLY_CODE, FISH_ID) %>% 
    dplyr::filter(n < nloci)
  
  # what loci are these indvs missing?
  missing_indvs_loci <- dataAll %>% 
    dplyr::distinct(SILLY_CODE, FISH_ID) %>% 
    tidyr::unite(col = "SillySource", c("SILLY_CODE", "FISH_ID"), sep = "_", remove = FALSE) %>% 
    tibble::add_column(!!!purrr::set_names(x = rep(NA_real_, nloci), nm = loci)) %>% 
    tidyr::pivot_longer(cols = c(-SILLY_CODE, -FISH_ID, -SillySource), names_to = "LOCUS", values_to = "na") %>% 
    dplyr::select(-na) %>% 
    dplyr::anti_join(dplyr::select(.data = dataAll, SILLY_CODE, FISH_ID, LOCUS), by = c("SILLY_CODE", "FISH_ID", "LOCUS")) %>% 
    tidyr::nest(missing_loci = LOCUS) %>% 
    dplyr::right_join(missing_indvs, by = c("SILLY_CODE", "FISH_ID")) %>% 
    dplyr::rename(n_loci_genotyped = n)
  
  names(missing_indvs_loci$missing_loci) <- missing_indvs_loci$SillySource
  
  # replace no calls (0's) with NA
  dataAll <- dataAll %>% 
    dplyr::mutate(ALLELE_1 = dplyr::na_if(ALLELE_1, "0"),
                  ALLELE_2 = dplyr::na_if(ALLELE_2, "0"))
  
  # filter out individuals missing loci, 
  dataAll_no_missing <- dataAll %>% 
    dplyr::anti_join(dplyr::select(.data = missing_indvs, SILLY_CODE, FISH_ID), by = c("SILLY_CODE", "FISH_ID")) 
  
  # what sillys have complete data?
  complete_sillys <- dataAll_no_missing$SILLY_CODE %>% 
    unique() %>% 
    sort()
  
  # What sillys had no individuals with complete data
  missing_sillys_indv <- setdiff(setdiff(sillyvec, missing_sillys), complete_sillys)

  dataAll_sillys <- dataAll$SILLY_CODE %>% 
    unique() %>% 
    sort()
  
  # did all indvs from a silly get dropped?
  missing_sillys_indv <- setdiff(setdiff(sillyvec, missing_sillys), dataAll_sillys)
  
  lapply(dataAll_sillys, function(silly){ 
    
    message0 <- paste0(silly, ".gcl created ", match(silly, dataAll_sillys)," of ", length(dataAll_sillys), " completed.") 
    
    sillydata <- dataAll %>% 
      dplyr::filter(SILLY_CODE==silly)
  
    silly_df_cols <- rep(NA_real_, nloci*2) %>% 
      purrr::set_names(c(rbind(loci, paste0(loci, ".1")))) 
    
    silly_df0 <- sillydata %>%
      dplyr::arrange(LOCUS) %>%
      tidyr::pivot_longer(cols = c("ALLELE_1", "ALLELE_2"), values_to = "Allele") %>% 
      dplyr::mutate(scores_header = dplyr::case_when(name == "ALLELE_2" ~ paste0(LOCUS, ".1"), 
                                                     TRUE ~ LOCUS)) %>% 
      dplyr::select(-LOCUS, -name) %>% 
      tidyr::pivot_wider(names_from = scores_header, values_from = Allele, names_sep = "" ) %>% 
      dplyr::mutate(
        CAPTURE_DATE = lubridate::as_date(CAPTURE_DATE),
        END_CAPTURE_DATE = lubridate::as_date(END_CAPTURE_DATE),
        SillySource = paste(SILLY_CODE, FISH_ID, sep = "_")
      )
    
    silly_df <- tibble::add_column(silly_df0, !!!silly_df_cols[setdiff(names(silly_df_cols), names(silly_df0))]) %>%
      dplyr::select(
        FK_FISH_ID = FISH_ID,
        COLLECTION_ID,
        SILLY_CODE,
        PLATE_ID,
        PK_TISSUE_TYPE,
        CAPTURE_LOCATION,
        CAPTURE_DATE,
        END_CAPTURE_DATE,
        MESH_SIZE,
        MESH_SIZE_COMMENT,
        LATITUDE,
        LONGITUDE,
        AGENCY,
        VIAL_BARCODE,
        DNA_TRAY_CODE,
        DNA_TRAY_WELL_CODE,
        DNA_TRAY_WELL_POS,
        CONTAINER_ARRAY_TYPE_ID,
        SillySource,
        tidyselect::all_of(names(silly_df_cols))
      ) %>%
      dplyr::arrange(FK_FISH_ID)
    
    message(message0)
    
    assign(paste0(silly, ".gcl"), silly_df, pos = 1, .GlobalEnv)
    
  })
  
  if(length(missing_sillys) >= 1){
    
    warning(paste0("The following sillys had no data in LOKI for any of the loci in LocusControl:\n", paste0(missing_sillys, collapse = "\n")), call. = FALSE)
    
  }
  
  if(length(missing_sillys_indv) >= 1){
    
    warning(paste0("The following sillys had no individuals with complete data in LOKI for the loci in LocusControl:\n", paste0(missing_sillys_indv, collapse = "\n")), call. = FALSE)
    
  }
  
  n_missing <- missing_indvs_loci %>% 
    dplyr::mutate(n_loci_missing = nloci - n_loci_genotyped) %>% 
    dplyr::count(SILLY_CODE, n_loci_missing) %>%
    dplyr::rename(n_indv = n) %>% 
    dplyr::mutate(silly_n_miss = paste0(SILLY_CODE, " (", n_indv, " individuals missing ", n_loci_missing, " loci)")) %>% 
    dplyr::pull(silly_n_miss)
  
  if(nrow(missing_indvs_loci) >= 1){
    
    warning(paste0("The following sillys had individuals with missing data for one or more loci:\n", paste(n_missing, collapse = "\n"), "\n"), call. = FALSE)
    
    warning(paste0("A table of loci missing data for each individual has been assigned to the object 'missing_indvs_loci'\n"), call. = FALSE)
    
    assign(x = "missing_indvs_loci", value = missing_indvs_loci, pos = 1, .GlobalEnv)
    
  } 
  
  if(!nrow(missing_indvs_loci) >= 1){
    
    print("The *.gcl objects created have data for all loci in LocusControl\n")
    
  }
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime)
  
}