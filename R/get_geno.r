#' Get Genotypes
#'
#' This function pulls a genotypes report from LOKI and writes the data to an Excel file.
#'
#' @param project_name A character vector of GCL project names to pull from LOKI (default = NULL).
#' @param sillyvec A character vector of silly codes without the ".gcl" extension (default = NULL).
#' @param loci A character vector of locus names as spelled in LOKI (default = NULL).
#' @param file The full file path with .xlsx extension for writing output genotypes file; with "\\" or "/" separator between folder.
#' @param username Your user name for accessing LOKI.
#' @param password Your password for accessing LOKI.
#' @param open.file Logical, if set to TRUE, the .csv file will open after it has been written (default = FALSE).
#' @param proj.stats whether you want a formatted Excel workbook with project genotype statistics (TRUE) or an unformatted Excel file without statistics (FALSE) (default = FALSE)
#' 
#' @returns A slim, tibble and a UTF-8 coded .csv file with single column genotypes in a wide format (1 column per locus):
#'     \itemize{
#'       \item \code{LAB_PROJECT_NAME}: the lab project name
#'       \item \code{FK_COLLECTION_ID}: the collection ID (internal to LOKI)
#'       \item \code{SILLY_CODE}: the silly code
#'       \item \code{FK_FISH_ID}: the fish ID
#'       \item \code{DNA_PLATE_ID}: the plate ID that was genotyped (only if `project_name` is given) 
#'       \item \code{Locus_1}: genotype for the first locus, 1 column (i.e., all alleles, separated by "/")
#'       \item \code{Locus_n}: genotype for the nth locus, 1 column (i.e., all alleles, separated by "/")
#'       }
#' 
#' @details This function requires the R package "RJDBC" for connecting to LOKI.
#' It connects to LOKI using the provided \code{username} and \code{password}, and retrieves genotypes and and slim attributes.
#' The genotypes report will only contain DNA plate IDs if project_name is supplied. When `project_name` is supplied and `proj.stats = TRUE`, the function calls on `GCLr::proj_geno_stats()` to produce a formatted Excel workbook with genotype statistics.
#' 
#' @seealso [GCLr::proj_geno_stats()]
#'
#' @examples
#' \dontrun{
#' Geno_locisilly <- GCLr::get_geno(sillyvec = sillyvec, loci = loci, file = path.expand("~/TestGenotypesReport_sillyloci.xlsx"), username = "awbarclay", password = password, project_name = NULL)
#' Geno_silly <- GCLr::get_geno(sillyvec = sillyvec, loci = NULL, file = path.expand("~/TestGenotypesReport_silly.xlsx"), username = "awbarclay", password = password, project_name = NULL)
#' Geno_proj <- GCLr::get_geno(project_name = "K205", file = path.expand("~/TestGenotypesReport_proj.xlsx"), username = "awbarclay", password = password, open.file = TRUE)
#' }
#' 
#' @export
get_geno <- function(project_name = NULL, sillyvec = NULL, loci = NULL, file, username, password, open.file = FALSE, proj.stats = FALSE){  

  # Recording function start time
  start.time <- Sys.time()
  
  #Make sure a project name is supplied if proj.stats = TRUE
  if(proj.stats == TRUE & is.null(project_name)){
    
    stop("A project name must be supplied when proj.stats = TRUE")
    
  }
  
  # Checking to make sure a file path is supplied
  if(!exists("file") | !length(grep("*.xlsx", file, value = FALSE)) == 1){
    
    stop("The user must supply a file path with xlsx extension for writing out genotypes table.")
    
  }
  
  # Checking to make sure the correct combination of arguments is being use. If wrong, the function will stop and print an error message to the console.
  if(is.null(sillyvec) & is.null(loci) & is.null(project_name)|
     !is.null(sillyvec) & !is.null(loci) & !is.null(project_name)|
     !is.null(sillyvec) & is.null(loci) & !is.null(project_name)|
     is.null(sillyvec) & !is.null(loci) & !is.null(project_name)){ 
    
    stop("The user must supply one of the following argument combinations: sillyvec (for all loci and individuals for each silly), 
         sillyvec and loci (all individuals for supplied locus list), or 
         project_name (for all individuals and loci in a given project)")
    
    }
  
  
  # Setting default java.parameters
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
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
    
    dataAll <- dplyr::left_join(geno, proj, by = c("LAB_PROJECT_NAME", "SILLY_CODE", "FK_FISH_ID" = "FISH_NO"))
    
    dataAll <- dataAll %>% 
      dplyr::mutate(DNA_PLATE_ID = PLATE_ID) %>% 
      dplyr::filter(LAB_PROJECT_NAME %in% project_name) %>% 
      dplyr::select(-ALLELE_1, -ALLELE_2) %>%
      dplyr::rename(ALLELE_1 = ALLELE_1_FIXED, ALLELE_2 = ALLELE_2_FIXED) %>% 
      tidyr::unite(GENO, ALLELE_1, ALLELE_2, sep = "/", remove = FALSE) %>% 
      dplyr::mutate(ALLELES = dplyr::case_when(PLOIDY == "D" ~ GENO,
                                               PLOIDY == "H" ~ ALLELE_1)) %>% 
      dplyr::select(LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, DNA_PLATE_ID, LOCUS, ALLELES) %>%
      tidyr::pivot_wider(names_from = LOCUS, values_from = ALLELES)
     
  } else {
    
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
  discon <- RJDBC::dbDisconnect(con)
  
  # Checking to see if no data was pulled for a given project. If true, the function will stop and print an error message in the console.
  if(nrow(dataAll) == 0 & exists("project_name")){
    
    stop("No data exist for the supplied lab project name. Check to make sure the name spelled correctly. Also, project names for genotypes imported prior to October 2016 may not be included in LOKI lookup table.")
   
     }
    
  if(proj.stats == TRUE){
    
    GCLr::proj_geno_stats(report = dataAll, file = file)
    
  }else{
    
    # Write data to an excel file
    dataAll %>% 
      writexl::write_xlsx(path = file)
    
  }
    
  # Open file in excel if open.file is set to TRUE
  if(open.file == TRUE){
    
    shell(file, wait = FALSE)
    
    }
  
  # Calculate the time it took the function to pull and write data and print time to console
  stop.time <- Sys.time() 
  
  message(paste0("Excel file writen to:", tools::file_path_as_absolute(file)))
  print(stop.time-start.time)
  
  # Return data to assign to an object
  return(dataAll)
   
}