#' Select High Quality QC Samples
#' 
#' This function pulls project genotypes from Loki and removes samples with genotypes at less than 80% of loci, then selects QC samples by plate (see details).
#' 
#' @param project the project name as it is spelled in Loki.
#' 
#' @param username Your user name for accessing LOKI.
#' 
#' @param password Your password for accessing LOKI.
#' 
#' @param ncores A numeric value for the number of cores to use in a \pkg{foreach} `%dopar%` loop (default = 4). 
#'
#' @details
#' When selecting samples for QC analysis, the function tries to select the normal QC card positions first. (i.e., `c("A3", "B4", "C5", "D6", "E7", "F8", "G9", "H10")`)
#' If any of those samples were removed for missing genotypes, the function selects a random sample from the same row as a replacement.  No samples are selected from extraction plates with less than 3 columns. 
#' 
#' @returns a tibble of all project samples containing the following variables: `Container Type`, `WBID`, `Barcode`, `Well Position`, `SILLY`, `Fish ID`, `Tissue`, `SILLY_Fish_ID`, `select`
#'          The function also writes out a csv file of the selection tibble (e.g., K215_QC_selection.csv) and another file of samples that were not removed for missing genotypes (e.g., K215_good_samples.csv) in case manual adjustments need to be made to the selection.
#'
#' @examples
#' \dontrun{
#' 
#'  select_qc_samps(project = "K215", username = "awbarclay", password = "mypassword", ncores = parallel::detectCores())
#'  
#' }
#'
#' @export
select_qc_samps <- function(project, username, password, ncores = 4){
  
  if(exists("LocusControl")){
    
    LocusCtl <- LocusControl
    
    rm(LocusControl, envir = .GlobalEnv)
    
    }
  
  suppressMessages(GCLr::loki2r_proj(project_name = project, username = username, password = password)) #Get project data
  
  #Put all samples in one object before removing bad samples.
  proj_samps0 <- lapply(project_sillys, function(silly){
    
    get(paste0(silly, ".gcl"))
    
  }) %>% dplyr::bind_rows()
  
  removed <- suppressMessages(GCLr::remove_ind_miss_loci(sillyvec = project_sillys)) # Remove samples with genotypes at < 80% of loci
  
  #Put all samples in one object after removing bad samples
  good_samps0 <- lapply(project_sillys, function(silly){
    
    get(paste0(silly, ".gcl"))
    
  }) %>% dplyr::bind_rows()
  
  #Get the plate numbers
  plates <- good_samps0$PLATE_ID %>% 
    unique()
  
  ### Get container types ###
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  qry <- paste0("SELECT * FROM AKFINADM.LU_CONTAINER_ARRAY_TYPE")

  container_key <- RJDBC::dbGetQuery(con, qry)  #Pulling data from LOKI using the connection and tissue data query
 
  discon <- RJDBC::dbDisconnect(con) # Disconnect from Loki
  
  ### End get container types ###

  tissue_info <- GCLr::get_tissue_data(sillyvec = project_sillys, username = username, password = password, import.vars = FALSE) %>% 
    tidyr::unite(col = "SillySource", SILLY_CODE, FK_FISH_ID, remove = FALSE)#Getting tissue data for workbench IDs
  
  extraction_info <- GCLr::get_extraction_info(plate_ids = plates, username = username, password = password) # Getting extraction info for well numbers.
  
  Well_pos_key <- tibble::tibble(Well_Code = sapply(LETTERS[1:8], function(letter){paste0(letter, 1:6)}) %>% 
    t() %>% 
    as.vector(), Well_Position = 1:48)
  
  good_samps <- dplyr::left_join(good_samps0 %>% dplyr::select(FK_FISH_ID, COLLECTION_ID, SILLY_CODE, PLATE_ID, PK_TISSUE_TYPE, SillySource), extraction_info, by = c("SILLY_CODE", "FK_FISH_ID" = "FISH_NO", "PLATE_ID" = "FK_PLATE_ID")) %>% 
    dplyr::left_join(tissue_info, by = c("SILLY_CODE", "FK_FISH_ID", "SillySource")) %>% 
    dplyr::select(Extraction_Well_Pos = WELL_NO, PLATE_ID, Container_Type = CONTAINER_ARRAY_TYPE_ID, WBID = DNA_TRAY_WORKBENCH_ID, Barcode = DNA_TRAY_CODE, Well_Code = DNA_TRAY_WELL_POS,  SILLY = SILLY_CODE, Fish_ID = FK_FISH_ID, Tissue = TISSUETYPE, SILLY_Fish_ID = SillySource) %>% 
    dplyr::left_join(Well_pos_key, by = "Well_Code") %>% 
    dplyr::left_join((tibble::as_tibble(container_key) %>% dplyr::select(LU_CONTAINER_ARRAY_TYPE_ID, CONTAINER_ARRAY_TYPE)), by = c("Container_Type" = "LU_CONTAINER_ARRAY_TYPE_ID")) %>% 
    dplyr::arrange(Fish_ID)
    
  proj_samps <- dplyr::left_join(proj_samps0 %>% dplyr::select(FK_FISH_ID, COLLECTION_ID, SILLY_CODE, PLATE_ID, PK_TISSUE_TYPE, SillySource), extraction_info, by = c("SILLY_CODE", "FK_FISH_ID" = "FISH_NO", "PLATE_ID" = "FK_PLATE_ID")) %>% 
    dplyr::right_join(tissue_info, by = c("SILLY_CODE", "FK_FISH_ID", "TISSUETYPE" = "PK_TISSUE_TYPE", "SillySource")) %>% 
    dplyr::select(Extraction_Well_Pos = WELL_NO, PLATE_ID, Container_Type = CONTAINER_ARRAY_TYPE_ID, WBID = DNA_TRAY_WORKBENCH_ID, Barcode = DNA_TRAY_CODE, Well_Code = DNA_TRAY_WELL_POS,  SILLY = SILLY_CODE, Fish_ID = FK_FISH_ID, Tissue = TISSUETYPE, SILLY_Fish_ID = SillySource) %>% 
    dplyr::left_join(Well_pos_key, by = "Well_Code") %>% 
    dplyr::left_join((tibble::as_tibble(container_key) %>% dplyr::select(LU_CONTAINER_ARRAY_TYPE_ID, CONTAINER_ARRAY_TYPE)), by = c("Container_Type" = "LU_CONTAINER_ARRAY_TYPE_ID")) %>% 
    dplyr::arrange(Fish_ID)
  
  #Select samples
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  # Start parallel loop
  
  `%dopar%` <- foreach::`%dopar%`
  
  selection <- foreach::foreach(plate = plates, .packages = c("tidyverse")) %dopar% {
    
    plate_samps <- good_samps %>% 
      dplyr::filter(PLATE_ID == plate) 
  
    ncols <- (as.numeric(substring(plate_samps$Extraction_Well_Pos, 2, nchar(plate_samps$Extraction_Well_Pos))) %>% 
               unique())[-c(1,2, 11, 12)] %>% length()
    
    select_well_codes <- c("A3", "B4", "C5", "D6", "E7", "F8", "G9", "H10") %>% 
      set_names(LETTERS[1:8])
    
    select0 <- plate_samps %>% 
      dplyr::filter(Extraction_Well_Pos %in% select_well_codes)
    
    if(dim(select0)[1] < ncols){
      
      missing_wells <- setdiff(select_well_codes[1:ncols], select0$Extraction_Well_Pos )
      
      missing_rows <- stringr::str_extract(missing_wells, pattern = "[:alpha:]")
    
      select_1 <- lapply(missing_rows, function(row){
        
        plate_samps %>% 
          dplyr::filter(grepl(pattern = row, x = Extraction_Well_Pos)) %>% 
                dplyr::sample_n(size = 1)
        
      }) %>% dplyr::bind_rows()
      
      select <- dplyr::bind_rows(select0, select_1)
      
    }else{select <- select0}
    
     select 
      
    } %>% 
    dplyr::bind_rows() %>% 
    dplyr::arrange(Fish_ID) %>% 
    dplyr::select(SILLY_Fish_ID) %>% 
    dplyr::mutate(select = 1)
  
  parallel::stopCluster(cl)  # End parallel loop
  
  # Write out all samples with selection column
  section_file <- paste0(project, "_QC_selection.csv")
  
  message(paste0("Writing out a file of all project silly samples with the selected samples indicated by the select column: ", getwd(), "/", section_file))
  
  select_out <- dplyr::left_join(proj_samps, selection, by = "SILLY_Fish_ID") %>%
    dplyr::select(`Container Type` = CONTAINER_ARRAY_TYPE, WBID, Barcode, `Well Position` = Well_Position, SILLY, `Fish ID` = Fish_ID, Tissue, SILLY_Fish_ID, select) %>% 
    dplyr::arrange(SILLY, `Fish ID`)
    
  readr::write_csv(select_out, section_file, na = "")
  
  # Write out only the samples that have sufficient genotypes
  good_samps_file <- paste0(project, "_good_samples.csv")
  
  message(paste0("Writing out a file containing samples that had genotypes at > 80% of loci to the following file: ", getwd(), "/", good_samps_file))
  
  goods_samples <- good_samps %>% 
    dplyr::select(`Container Type` = CONTAINER_ARRAY_TYPE, `WBID`, Barcode, `Well Position` = Well_Position, SILLY, `Fish ID` = Fish_ID, Tissue, SILLY_Fish_ID) %>% 
    dplyr::arrange(SILLY, `Fish ID`)
  
  readr::write_csv(goods_samples, good_samps_file, na = "")
  
  if(exists("LocusCtl")){
    
    assign(x = "LocusControl", value = LocusCtl, envir = .GlobalEnv)
    
    }
  
  return(select_out)
  
}