#' @title Read Project Genotypes
#'
#' @description This function is intended for use in the qc.R script. It pulls project genotypes to create "slim" .gcl objects for each silly and the LocusControl for all loci used in the project.
#' 
#' **Warning**: Genotypes (pre-October 2016) may not exist in the LOKI lookup table with a project name, so genotypes will have to be pulled by sillyvec and loci.
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
#' This function pulls project genotypes from LOKI database to create "slim" .gcl objects for each silly,
#' along with the LocusControl for all loci used in the project.
#' The function checks the combination of arguments to ensure the correct usage.
#' It requires the `tidyverse`, `RJDBC`, and `tools` packages.
#' The function also requires the `ojdbc8.jar` file, which will be copied to the appropriate location if it doesn't exist.
#'
#' @aliases read_project_genotypes.GCL.R, read_proj_geno
#' @export
#' 
#' @examples
#' read_proj_geno(project_name = c("P014", "P015", "P016"), sillyvec = c("SCIMA18", "SCIMA17"), 
#'                loci = c("One_E2", "One_MHC2_251", "One_Cytb_17"), username = "awbarclay", password = "password")
#'
#' @author Kyle Shedd, Andy Barclay
#' @updated 1/11/19 (Kyle Shedd), 10/5/18 (Andy Barclay), 4/15/19 (Andy Barclay)
read_proj_geno <- function(project_name = NULL, sillyvec = NULL, loci = NULL, username, password) {
  # Warning for use of deprecated function name.
  if (match.call()[[1]] != quote(read_proj_geno)) {
    warning("The function name 'read_project_genotypes.GCL.R' is deprecated. Please use 'read_proj_geno' instead.")
  }
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
  
  # Requiring the tidyverse and RJDBC packages. Packages will be installed if the user doesn't have them installed
  while (!require(tidyverse)) {install.packages("tidyverse") }
  while (!require(RJDBC)) {install.packages("RJDBC") }
  while (!require(tools)) {install.packages("tools") }
  
  # Making sure the ojdbc jar file exists on the users C drive. 
  # If no jar file exists it will be copied from the v drive to the appropriate location on the users C drive.
  if (!file.exists(path.expand("~/R"))) {
    dir <- path.expand("~/R")
    dir.create(dir)
    bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar",to = path.expand("~/R/ojdbc8.jar"))
  } else {
    if (!file.exists(path.expand("~/R/ojdbc8.jar"))) {
      bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar",to = path.expand("~/R/ojdbc8.jar"))
    }
  }
  
  # Setting default java.parameters
  options(java.parameters = "-Xmx10g")
  
  if (file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc8.jar", " ")  # https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
  } else {
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = path.expand("~/R/ojdbc8.jar"), " ")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Connect to LOKI
  
  url <- loki_url()
  
  con <- dbConnect(drv,url = url,user = username,password = password)
  
  #~~~~~~~~~~~~~~~~
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
      filter(LAB_PROJECT_NAME == project_name)
  }
  
  # Get list of unique sillys and assign `ProjectSillys` this is needed for qc script
  sillyvec <- unique(dataAll$SILLY_CODE)
  assign(x = "ProjectSillys", value = sillyvec, pos = 1)
  
  # Disconnect from LOKI
  discon <- RJDBC::dbDisconnect(con)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get list of unique loci and assign `LocusControl` this is needed for qc script
  loci <- sort(unique(dataAll$LOCUS))
  nloci <- length(loci)
  
  create_locuscontrol(locusnames = loci, username = username, password = password)
  
  assign(x = "loci", value = LocusControl$locusnames, pos = 1)
  assign(x = "nalleles", value = LocusControl$nalleles, pos = 1)
  assign(x = "ploidy", value = LocusControl$ploidy, pos = 1)
  assign(x = "alleles", value = LocusControl$alleles, pos = 1)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Join genotypes with extraction information (PLATE_ID) and make .gcl objects
  # This still needs work because it pulls all plates per fish, not just the project plates!!!
  
  attnames <- c("FK_FISH_ID", "COLLECTION_ID", "SILLY_CODE", "PLATE_ID", "PK_TISSUE_TYPE", "SillySource")
  
  data_master <- dataAll %>% 
    tidyr::unite(SillySource, SILLY_CODE, FK_FISH_ID, sep = "_", remove = FALSE) %>% 
    dplyr::select(-ALLELE_1, -ALLELE_2) %>%
    dplyr::rename(ALLELE_1 = ALLELE_1_FIXED, ALLELE_2 = ALLELE_2_FIXED) %>% 
    tidyr::unite(GENO, ALLELE_1, ALLELE_2, sep = "/", remove = FALSE) %>% 
    dplyr::mutate(ALLELES = dplyr::case_when(PLOIDY == "D" ~ GENO,
                                             PLOIDY == "H" ~ ALLELE_1)) %>% 
    dplyr::select(c(dplyr::all_of(attnames), "LOCUS", "ALLELE_1", "ALLELE_2", "ALLELES")) #%>% 
  # dplyr::select(-PLOIDY) %>% 
  # tidyr::spread(key = LOCUS, value = ALLELES)
  
  message("Data successfully pulled from LOKI, building SILLY.gcl objects\n")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Make .gcl objects by silly
  for (silly in sillyvec) {  
    
    data_silly <- dplyr::filter(data_master, SILLY_CODE == silly)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Adaptation of LOKI2R to allow pulling all SILLYs just based on ProjectID
    if (nrow(data_silly) == 0) {
      warning(paste(silly, "was not found in LOKI, it may not have been imported!!!", sep = " "))
      message(paste0(silly,".gcl NOT FOUND IN LOKI ", match(silly,sillyvec)," of ",length(sillyvec)," FAILED!!!"))
      next
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ids <- as.character(sort(unique(data_silly$FK_FISH_ID)))
    nind <- length(ids)
    
    if (!nind) {
      message0 <- paste0(silly," is empty.")
      message(message0)
      next
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Build scores array
    
    scores <- array(data = NA,
                    dim = c(nind, nloci, max(ploidy)), 
                    dimnames = list(ids, loci, paste0("Dose", seq(max(ploidy))))
    )
    
    sillyloci <- sort(unique(data_silly$LOCUS))  # These are the loci available for the current silly, this is needed to subset the tapply
    
    scores_allele_1 <- data_silly %>% 
      dplyr::filter(LOCUS %in% sillyloci) %>% 
      dplyr::select(FK_FISH_ID, LOCUS, ALLELE_1) %>% 
      dplyr::mutate(FK_FISH_ID = as.character(FK_FISH_ID)) %>% 
      tidyr::spread(LOCUS, ALLELE_1) %>% 
      tibble::column_to_rownames(var = "FK_FISH_ID") %>% 
      as.matrix()
    
    scores[ids, sillyloci, "Dose1"] <- scores_allele_1[ids, sillyloci]
    
    scores_allele_2 <- data_silly %>% 
      dplyr::filter(LOCUS %in% sillyloci) %>% 
      dplyr::select(FK_FISH_ID, LOCUS, ALLELE_2) %>% 
      dplyr::mutate(FK_FISH_ID = as.character(FK_FISH_ID)) %>% 
      tidyr::spread(LOCUS, ALLELE_2) %>% 
      tibble::column_to_rownames(var = "FK_FISH_ID") %>% 
      as.matrix()
    
    scores[ids, sillyloci, "Dose2"] <- scores_allele_2[ids, sillyloci]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Build counts array
    
    counts <- array(data = NA, 
                    dim = c(nind, nloci, max(nalleles)), 
                    dimnames = list(ids, loci, paste0("Allele ", seq(max(nalleles))))
    )
    
    for (ind in ids) {
      for (locus in loci) {
        for (al in seq(nalleles[locus])) {
          counts[ind, locus, al] <- sum(scores[ind, locus, seq(ploidy[locus])] == alleles[[locus]][al])
        }  # al
      }  # locus
      counts0 = counts[ind, , ]
      counts0[is.na(counts0)] <- 0
      zeroBOOL <- apply(counts0,1,sum) != ploidy
      counts[ind, zeroBOOL, ] <- NA
    }  # ind
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Build attributes data.frame
    
    attributes <- data_silly %>% 
      dplyr::select(dplyr::all_of(attnames)) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(FISH_ID = as.character(FK_FISH_ID),
                    PLATE_ID = as.character(PLATE_ID)) %>% 
      dplyr::arrange(FK_FISH_ID) %>% 
      tibble::column_to_rownames(var = "FISH_ID")
    
    message(paste0(silly, ".gcl created ", match(silly, sillyvec), " of ", length(sillyvec), " completed."))
    
    assign(paste0(silly, ".gcl"), list(counts = counts, scores = scores, n = nind, attributes = attributes), pos = 1)
    
  }  # silly
  
  stop.time <- Sys.time()
  fulltime <- stop.time - start.time
  print(fulltime) 
}
