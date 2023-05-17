#' This function connects to the Gene Conservation Lab Oracle database Loki and creates a "gcl" object containing genotypes and paired data for each sample.
#'
#' An object, similar to the rubias input format, is created for each collection (a.k.a silly) in sillyvec containing genotypes for each locus in LocusControl$locusnames. The default for this function is to only include fish that have been genotyped for all loci in LocusControl$locusnames; however, the function will include all fish when include_missing = TRUE, which is not recommended or needed for most analyses.
#'
#' @param sillyvec a character vector of silly codes (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART09"))
#' @param username your Oracle user name
#' @param password your Oracle password (contact Eric Lardizabal if you don't have a username or password)
#' @param test_type the test type ("SNP", "GTSNP", "MSAT) you would like to pull from Loki - default is test_type = "SNP". This argument is a workaround. Some collections in Loki contain data for both SNP and GTSNP test types. If genotypes for both test types were pulled at the same time for the same individuals and loci, this would cause the function to stop and throw an error message. If you want data for both test types in your workspace, you will need to run this function separately for each test type. If you pull data for both test types for the same silly code, you will need to rename your '.gcl' object for the first test type before pulling data for the second test type, otherwise the first '.gcl' object will be overwritten.
#' @param include_missing logical statement whether to include all samples even if they were never genotyped for all loci in LocusControl.
#'
#' @export
#'
#' @return   This function returns a tibble of two-column genetic data with 19 columns of sample attributes preceding it.
#'
#' @note This function requires a LocusControl object. Run CreateLocusControl.GCL prior to this function.
#'
#' @examples
#'
#' CreateLocusControl.GCL(markersuite = "Sockeye2011_96SNPs", username ="awbarclay", password = "mypassword")
#'
#' sillyvec <- c("SUCIWS06", "SUCIWS07", "SUCIWS08", "SUCIWS09", "SUCIWS10", "SUCIWS11", "SUCIWS12", "SUCIWS13", "SCIMA22")
#'
#' LOKI2R.GCL(sillyvec = sillyvec, username = "awbarclay", password = .password, test_type = "SNP", include_missing = TRUE)
#'
LOKI2R.GCL <- function(sillyvec, username, password, test_type = c("SNP", "GTSNP", "MSAT")[1], include_missing = FALSE){

  if(!exists("LocusControl")){

    stop("'LocusControl' not yet built.")

  }

  if(length(test_type)>1){

    stop("Only one test_type can be supplied at a time.")

  }

  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(RJDBC, tidyverse, lubridate) #Install packages, if not in library and then load them.

  # This copies the "odbc8.jar" file to the R folder on your computer if it doesn't exist there. This file contains the java odbc drivers needed for RJDBC

  if(!file.exists(path.expand("~/R"))){

    dir <- path.expand("~/R")

    dir.create(dir)

    bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))

  } else {

    if(!file.exists(path.expand("~/R/ojdbc8.jar"))){

      bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))

    }

  }

  start.time <- Sys.time()

  options(java.parameters = "-Xmx10g")

  if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {

    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc8.jar", " ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib

  } else {

    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = path.expand("~/R/ojdbc8.jar"), " ")

  }

  url <- LOKI_URL.GCL() #This is a function that gets the correct URL to access the database on the oracle cloud

  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password) #The database connection

  loci <- LocusControl$locusnames

  nloci <- length(loci)

  ploidy <- LocusControl$ploidy

  hap_loci <- ploidy[ploidy==1] %>% names()

  alleles <- LocusControl$alleles

  nalleles <- LocusControl$nalleles

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
      filter(TEST_TYPE == "GTSNP") %>%
      pull(SILLY_CODE) %>%
      unique()

  }

  if(test_type == "GTSNP" & sum(test_types == "SNP")>0){

    SNP_sillys <- dataAll %>%
      filter(TEST_TYPE == "SNP") %>%
      pull(SILLY_CODE) %>%
      unique()

  }

  dataAll <- filter(dataAll, TEST_TYPE == test_type) %>%
    select(-TEST_TYPE)# Filter dataAll for the supplied test_type

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
      purrr::set_names(c(loci, paste0(loci, ".1")) %>% sort())

    silly_df0 <- sillydata %>%
      dplyr::arrange(LOCUS) %>%
      tidyr::pivot_longer(cols = c("ALLELE_1", "ALLELE_2"), values_to = "Allele") %>%
      dplyr::mutate(scores_header = case_when(name == "ALLELE_2" ~ paste0(LOCUS, ".1"),
                                              TRUE ~ LOCUS)) %>%
      dplyr::select(-LOCUS, -name) %>%
      tidyr::pivot_wider(names_from = scores_header, values_from = Allele, names_sep="" ) %>%
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
