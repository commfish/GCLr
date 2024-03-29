#' Create GCL Objects
#'
#' This function connects to LOKI and creates a .gcl object for each silly in `sillyvec` containing genotypes for each locus in LocusControl$locusnames.The default for this function is to only include fish that have been genotyped for all loci in `LocusControl$locusnames`; however, the function will include all fish when include_missing = TRUE, which is not recommended or needed for most analyses. If a silly has no fish with genotypes for the loci in `LocusControl`, no object will be created and a message will appear listing the sillys that had no data.
#'
#' @param sillyvec a character vector of silly codes without the .gcl extension.
#' 
#' @param username your state user name
#' 
#' @param password your password used to access LOKI; see Eric Lardizabal if you need to set up a password
#' 
#' @param test_type the test type ("SNP" or "GTSNP") you would like to pull from Loki. (default = "SNP")
#' 
#' @param include_missing whether to include all fish even if they were never genotyped for all loci in LocusControl (default = FALSE)
#' 
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)
#'
#' @return This function assigns a tibble with the following columns for each silly:
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
#' @note This function requires a LocusControl object. Run [GCLr::create_locuscontrol()] prior to this function.
#'    
#' @examples
#' \dontrun{
#'   sillyvec <- c("SUCIWS06", "SUCIWS07", "SUCIWS08", "SUCIWS09", "SUCIWS10", "SUCIWS11", "SUCIWS12", "SUCIWS13", "SCIMA22")
#'   loki2r(sillyvec = sillyvec, username = "awbarclay", password = .password, test_type = "SNP", include_missing = TRUE, LocusCtl = LocusControl)
#' }
#' 
#' @export            
loki2r <- function(sillyvec, username, password, test_type = c("SNP", "GTSNP", "MSAT")[1], include_missing = FALSE, LocusCtl = LocusControl){
  
  if(length(test_type) > 1){
    
    stop("Only one test_type can be supplied at a time.")
    
  }
  
  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  loci <- LocusCtl$locusnames
  
  nloci <- length(loci)
  
  ploidy <- LocusCtl$ploidy
  
  hap_loci <- ploidy[ploidy==1] %>% names()
  
  alleles <- LocusCtl$alleles
  
  nalleles <- LocusCtl$nalleles 
  
  gnoqry <- paste("SELECT * FROM AKFINADM.V_GNOQRY WHERE LOCUS IN (", paste0("'", loci, "'", collapse = ","), ") AND SILLY_CODE IN (", paste0("'", sillyvec, "'", collapse = ","), ")", sep = "") #Gentoype query
  
  dataAll0 <- RJDBC::dbGetQuery(con, gnoqry)  #Pulling data from LOKI using the connection and genotype query
  
  dataAllbool <- dataAll0$PK_TISSUE_TYPE == dataAll0$PREFERRED_TISSUE 
  
  bothNAbool <- is.na(dataAll0$PREFERRED_TISSUE) &  is.na(dataAll0$PK_TISSUE_TYPE)
  
  dataAllbool[is.na(dataAllbool)] <- FALSE 
  
  dataAllbool[bothNAbool] <- TRUE  
  
  dataAll <- dataAll0[dataAllbool, ] %>% 
    tibble::as_tibble() %>% 
    dplyr::select(-PREFERRED_TISSUE)
  
  discon <- RJDBC::dbDisconnect(con)
  
  # Are there any sillys that contain data for a TEST_TYPE other that the one supplied?
  # If so, a warning message will be produced listing out these sillys.
  test_types <- unique(dataAll$TEST_TYPE)
  
  if(sum(test_types %in% test_type) == 0){
    
    stop(paste0("The supplied sillys do not contain data for test_type = ", test_type))
    
  }
  
  if(test_type == "SNP" & sum(test_types == "GTSNP")>0){
    
    GTSNP_sillys <- dataAll %>% 
      dplyr::filter(TEST_TYPE == "GTSNP") %>% 
      dplyr::pull(SILLY_CODE) %>% 
      unique()
    
  }
  
  if(test_type == "GTSNP" & sum(test_types == "SNP")>0){
    
    SNP_sillys <- dataAll %>% 
      dplyr::filter(TEST_TYPE == "SNP") %>% 
      dplyr::pull(SILLY_CODE) %>% 
      unique()
    
  }
  
  dataAll <- dplyr::filter(dataAll, TEST_TYPE == test_type) %>% 
    dplyr::select(-TEST_TYPE)# Filter dataAll for the supplied test_type
    
  # what sillys have no data for any of these loci?
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
    dplyr::mutate(ALLELE_1 = gsub(pattern = "^0$", replacement = NA_character_, x = ALLELE_1),
                  ALLELE_2 = gsub(pattern = "^0$", replacement = NA_character_, x = ALLELE_2))
  
  # filter out individuals missing loci, 
  dataAll_no_missing <- dataAll %>% 
    dplyr::anti_join(dplyr::select(.data = missing_indvs, SILLY_CODE, FISH_ID), by = c("SILLY_CODE", "FISH_ID")) 

  # what sillys have complete data?
  complete_sillys <- dataAll_no_missing$SILLY_CODE %>% 
    unique() %>% 
    sort()
  
  # What sillys had no individuals with complete data
  missing_sillys_indv <- setdiff(setdiff(sillyvec, missing_sillys), complete_sillys)
  
  if(include_missing == FALSE){
    
    dataAll <- dataAll_no_missing
    
  }
  
  dataAll_sillys <- dataAll$SILLY_CODE %>% 
    unique() %>% 
    sort()
    
  # did all indvs from a silly get dropped?
  missing_sillys_indv <- setdiff(setdiff(sillyvec, missing_sillys), dataAll_sillys)
  
  lapply(dataAll_sillys, function(silly){ 
    
    message0 <- paste0(silly, ".gcl created ", match(silly, dataAll_sillys)," of ", length(dataAll_sillys), " completed.") 
    
    sillydata <- dataAll %>% 
      dplyr::filter(SILLY_CODE==silly)
    
    ids <- sillydata$FISH_ID %>% 
      unique() %>% 
      sort() %>% 
      as.character()
    
    sillyvials <- paste(silly, ids, sep = "_")
    
    nind <- length(sillyvials)
    
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
    
    #Make sure all haploid loci have NA's for allele 2
    if(length(hap_loci) > 0){
      
      silly_df[ paste0(hap_loci, ".1")] <- NA_character_
       
    }
    
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
  
  if(nrow(missing_indvs_loci) >= 1 & include_missing == FALSE){
    
    warning(paste0("The following sillys had individuals that were removed due to missing data for one or more loci:\n", paste(n_missing, collapse = "\n"), "\n"), call. = FALSE)
    
    warning(paste0("A table of loci missing data for each individual has been assigned to the object 'missing_indvs_loci'\n"), call. = FALSE)
    
    assign(x = "missing_indvs_loci", value = missing_indvs_loci, pos = 1, .GlobalEnv)
    
  } 
  
  if(nrow(missing_indvs_loci) >= 1 & include_missing == TRUE){
    
    warning(paste0("The following sillys had individuals with missing data for one or more loci:\n", paste(n_missing, collapse = "\n"), "\n"), call. = FALSE)
    
    warning(paste0("A table of loci missing data for each individual has been assigned to the object 'missing_indvs_loci'\n"), call. = FALSE)
    
    assign(x = "missing_indvs_loci", value = missing_indvs_loci, pos = 1, .GlobalEnv)
    
  } 
  
  if(!nrow(missing_indvs_loci) >= 1){
    
    print("The *.gcl objects created have data for all loci in LocusControl\n")
    
  }
  
  if(exists("GTSNP_sillys")){
    
    warning(paste0("GTSNP test_type data was removed before creating '.GCL' objects for the following silly codes: ", paste0(GTSNP_sillys, collapse = ", ")))
    
  }
  
  if(exists("SNP_sillys")){
    
    warning(paste0("SNP test_type data was removed before creating '.GCL' objects for the following silly codes: ", paste0(SNP_sillys, collapse = ", ")))
    
  }
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime)
  
}