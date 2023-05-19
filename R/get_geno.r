get_geno <- function(project_name = NULL, sillyvec = NULL, loci = NULL, path, username, password, open.file = FALSE){  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function pulls a genotypes report from LOKI and writes the data to aUTF-8 coded csv file.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #  sillyvec - character vector of SILLYs - example: c("SCIMA18","SCIMA17")
  #
  #  loci - character vector of locus names as they are spelled in LOKI - example: c("One_E2","One_MHC2_251","One_Cytb_17")
  #
  #  path - file path with csv extension for writing output file with "\\" separator between folders example: "V:\\Analysis\\Staff\\Andy Barclay\\R\Genotypes report"
  #
  #  username - your user name for accessing LOKI through R example: "awbarclay"
  #
  #  password - your password for accessing LOKI through R example: "password"
  #
  #  project_name - a character vector of one or more project names as spelled in LOKI example: c("P014","P015","P016")
  #
  #  open.file - logical, if set to TRUE, the CSV file will open after it has been written. Default is FALSE.
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #  CSV file - Crosstab format of silly individuals by loci with both alleles in one column. 
  #
  #  Function returns a tibble of the genotypes
  # 
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # password <- "************"
  # username <- "awbarclay"
  # sillyvec <- c("SCHIG09s","SCHIG09j","SCHIG11s","SCHIGL11s","SCHIG11j")
  # loci <- c("One_E2","One_MHC2_251","One_Cytb_17")
  # project_name <- c("P014","P015","P016")
  #
  # Geno_locisilly <- get_geno(sillyvec = sillyvec, loci = loci, path = "TestGenotypesReport_sillyloci.csv", username = username, password = password, project_name = NULL)
  #
  # Geno_silly <- get_geno(sillyvec = sillyvec, loci = NULL, path = "TestGenotypesReport_silly.csv", username = username, password = password, project_name = NULL)
  #
  # Geno_proj <- get_geno(sillyvec = NULL, loci = NULL, path = "TestGenotypesReport_proj.csv", username = username, password = password, project_name = "S180", open.file = TRUE)
  #
  # Geno_P014_16 <- get_geno(sillyvec = NULL, loci = NULL, path = "TestGenotypesReport_P014_16.csv", username = username, password = password, project_name = project_name)
  # 
  # Geno_olddata <- get_geno(sillyvec = sillyvec, loci = NULL, path = "TestGenotypesReport_olddata.csv", username = username, password = password, project_name = NULL, open.file = TRUE)
  #
  # Geno_loci <- get_geno(sillyvec = NULL, loci = loci, path = "TestGenotypesReport_allsillys.csv", username = username, password = password, project_name = NULL, open.file = TRUE)
  # 
  # Notes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # The genotypes report will only contain DNA plate IDs if project_name is supplied.
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Recording function start time
  start.time <- Sys.time()
  
  # Checking to make sure the correct combination of arguments is being use. If wrong, the function will stop and print an error message to the console.
  if(is.null(sillyvec) & is.null(loci) & is.null(project_name)|
     !is.null(sillyvec) & !is.null(loci) & !is.null(project_name)|
     !is.null(sillyvec) & is.null(loci) & !is.null(project_name)|
     is.null(sillyvec) & !is.null(loci) & !is.null(project_name)){ 
    
    stop("The user must supply one of the following argument combinations: sillyvec (for all loci and individuals for each silly), 
         sillyvec and loci (all individuals for supplied locus list), or 
         project_name (for all individuals and loci in a given project)")
    
    }
  
  # Checking to make sure a file path is supplied
  if(!exists("path") | !length(grep("*.csv", path, value = FALSE)) == 1){
    
    stop("The user must supply a file path with csv extension for writing out genotypes table.")
    
    }
  
  # Requiring the tidyverse and RJDBC packages. Packages will be installed if the user doesn't have them installed
  
  # Making sure the ojdbc jar file exists on the users C drive. 
  # If no jar file exists it will be copied from the v drive to the appropriate location on the users C drive.
  if(!file.exists(path.expand("~/R"))){
    
    dir <- path.expand("~/R")
    
    dir.create(dir)
    
    bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
    
     } else {
    
           if(!file.exists(path.expand("~/R/ojdbc8.jar"))){
      
               bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
      
         }
  
       }
  
  # Setting default java.parameters
  options(java.parameters = "-Xmx10g")
  
  if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
    
        drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc8.jar"," ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
  
        } else {
      
          drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = path.expand("~/R/ojdbc8.jar"), " ")
  
        }
  
  url <- loki_url()
  
  # Connect to LOKI
  con <- DBI::dbConnect(drv, url = url, user = username, password = password)
         
  # Creating java query when sillyvec and loci are supplied.  
  if(!is.null(sillyvec) & !is.null(loci)){
    
    gnoqry <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1, ALLELE_2, ALLELE_1_FIXED, ALLELE_2_FIXED FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE LOCUS IN (", paste0("'", loci, "'", collapse = ","), ") AND SILLY_CODE IN", "(", paste0("'", sillyvec, "'", collapse = ","), ")")
  
    } 
  
  # Creating java query when only sillyvec is supplied
  if(!is.null(sillyvec) & is.null(loci)){
    
    gnoqry <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1, ALLELE_2, ALLELE_1_FIXED, ALLELE_2_FIXED FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE SILLY_CODE IN", "(", paste0("'", sillyvec, "'", collapse = ","), ")")
  
    }
  
  # Creating java query when only project_name is supplied.
  if(!is.null(project_name)){
    
    gnoqry <- paste0("SELECT * FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO GENO WHERE EXISTS (SELECT * FROM AKFINADM.V_LAB_PROJECT_WELL LPW WHERE LPW.LAB_PROJECT_NAME IN (", paste0("'", project_name, "'", collapse = ","), ") AND LPW.SILLY_CODE = GENO.SILLY_CODE AND LPW.FISH_NO = GENO.FK_FISH_ID)")
  
    projqry <- paste0("SELECT * FROM AKFINADM.V_LAB_PROJECT_WELL WHERE LAB_PROJECT_NAME IN (", paste0("'", project_name, "'", collapse = ","), ")")
    
  }
  
  # Creating java query when only loci is supplied.
  if(is.null(sillyvec) & !is.null(loci)){
    
    gnoqry <- paste0("SELECT LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, PLOIDY, ALLELE_1, ALLELE_2, ALLELE_1_FIXED, ALLELE_2_FIXED FROM AKFINADM.V_GEN_TEST_RESULTS_BOTHGENO WHERE LOCUS IN (", paste0("'", loci, "'", collapse = ","), ")")
    
  } 
  
  # Pull genotypes and concatenate alleles into one column with "/" separator
   
  if(!is.null(project_name)) {
    
    geno <- RJDBC::dbGetQuery(con, gnoqry) %>% 
      tibble::as_tibble() 
    
    proj <- RJDBC::dbGetQuery(con, projqry) %>% 
      tibble::as_tibble() 
    
    dataAll <- left_join(geno, proj, by = c("LAB_PROJECT_NAME", "SILLY_CODE", "FK_FISH_ID" = "FISH_NO"))
    
    dataAll <- dataAll %>% 
      dplyr::mutate(DNA_PLATE_ID = PLATE_ID) %>% 
      dplyr::filter(LAB_PROJECT_NAME == project_name) %>% 
      dplyr::select(-ALLELE_1, -ALLELE_2) %>%
      dplyr::rename(ALLELE_1 = ALLELE_1_FIXED, ALLELE_2 = ALLELE_2_FIXED) %>% 
      tidyr::unite(GENO, ALLELE_1, ALLELE_2, sep = "/", remove = FALSE) %>% 
      dplyr::mutate(ALLELES = dplyr::case_when(PLOIDY == "D" ~ GENO,
                                               PLOIDY == "H" ~ ALLELE_1)) %>% 
      dplyr::select(LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, DNA_PLATE_ID, LOCUS, ALLELES) %>%
      tidyr::pivot_wider(names_from = LOCUS, values_from = ALLELES)
     
  }else{
    
    dataAll <- RJDBC::dbGetQuery(con, gnoqry) %>% 
      tibble::as_tibble() 
      
    dataAll <- dataAll %>% 
      dplyr::select(-ALLELE_1, -ALLELE_2) %>%
      dplyr::rename(ALLELE_1 = ALLELE_1_FIXED, ALLELE_2 = ALLELE_2_FIXED) %>% 
      tidyr::unite(GENO, ALLELE_1, ALLELE_2, sep = "/", remove = FALSE) %>% 
      dplyr::mutate(ALLELES = dplyr::case_when(PLOIDY == "D" ~ GENO,
                                               PLOIDY == "H" ~ ALLELE_1)) %>% 
      dplyr::select(LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, LOCUS, ALLELES) %>%
      tidyr::pivot_wider(names_from = LOCUS, values_from = ALLELES)
    
    }
  
   
  # Disconnect from LOKI
  discon <- DBI::dbDisconnect(con)
  
  # Checking to see if no data was pulled for a given project. If true, the function will stop and print an error message in the console.
  if(nrow(dataAll) == 0 & exists("project_name")){
    
    stop("No data exist for the supplied lab project name. Check to make sure the name spelled correctly. Also, project names for genotypes imported prior to October 2016 may not be included in LOKI lookup table.")
   
     }
    
  # Write data to an excel csv file
  dataAll %>% 
    readr::write_excel_csv(file = path)
  
  # Open CSV file in excel if open.file is set to TRUE
  if(open.file == TRUE){
    
    shell(path, wait = FALSE)
    
    }
  
  # Calculate the time it took the function to pull and write data and print time to console
  stop.time <- Sys.time() 
  
  message(paste0("CSV file writen to:", tools::file_path_as_absolute(path)))
  print(stop.time-start.time)
  
  # Return data to assign to an object
  return(dataAll)
   
}

