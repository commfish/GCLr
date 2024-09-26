#' Get Quality Control Loki Bounce Summary and Write to Excel Workbook
#' 
#' This function produces an Excel file containing the project genotypes, the proportion of project fish with at least 80% of loci with scores, average success rate by plate, conflicts and success rates for individual samples, and conflicts and success rates by locus 
#' 
#' @param QC_directory the project QC directory
#' @param project the project name; must be spelled exactly the way it is in iStrategy.
#' @param username Your user name for accessing LOKI.
#' @param password Your password for accessing LOKI.
#' 
#' @details This function requires R package "RJDBC" for connecting to LOKI and calls on `GCLr::loki2r_proj()`, `GCLr::get_geno()`, 
#' and `GCLr::combine_conflicts()` to get the information it needs to create a QC bounce Excel Workbook. In order for this function to work,
#' there must be a folder named 'Conflict Reports' in the supplied `QC_directory` that contains Loki Concordance report(s) produced by the genotypes importer.
#' All concordance reports will be read in by the function as long as they include "Concordance" in their file names, so make sure that the folder only includes the reports you want summarized.
#' 
#' @seealso [GCLr::loki2r_proj()]
#' @seealso [GCLr::get_geno()]
#' @seealso [GCLr::combine_conflicts()]
#' 
#' @returns This function writes an Excel workbook with five sheets: 
#' \itemize{
#'   \item \code{proj_Genotypes_Table}: a formatted genotypes table
#'   \item \code{Calculations}: number of fish with success rates < 80%, total number of project fish, and Proportion of fish with success rates greater than 80%
#'   \item \code{Avg success by plate}: a table with *PlateID* and *Average Success Rate* columns
#'   \item \code{Conflicts by individual}: a table with *SillySource* (aka fish ID), *PlateID*, and *success_rate* columns and associated columns for each type of conflict in the concordance file(s) (e.g., Homo-Het, Het-Homo, Homo-Homo), and a column of *Total_conflicts*
#'   \item \code{Conflicts by locus}: a table with *locus*, *PlateID*, and *success_rate* columns and associated columns for each type of conflict in the concordance file(s) (e.g., Homo-Het, Het-Homo, Homo-Homo), and a column of *Total_conflicts*
#'   }
#' 
#' #' @examples
#' \dontrun{
#' 
#' QC_directory <- "V:/Lab/Genotyping/SNP Projects/Sockeye/Project S266 Port Moller Inseason 2024/QC"
#' project <- "S266"
#' username <- "awbarclay"
#' password <- scan("~/R/usr.pw", what = "")[2]
#' 
#' qc_bounce(QC_directory, project, username, password)
#' 
#'  } 
#' @export
qc_bounce <- function(QC_directory, project, username, password){
  
  setwd(QC_directory)
  
  out.file <- paste0(QC_directory, "/", project, "_genotypes_report.xlsx")
  
  GCLr::loki2r_proj(project_name = project, username = username, password = password)
  
  genotypes <- GCLr::get_geno(project_name = project, sillyvec = NULL, loci = NULL, file = out.file, username, password, proj.stats = TRUE)
  
  #Get sillyvec
  sillyvec <- genotypes$SILLY_CODE %>% 
    unique()
  
  loci <- names(genotypes)[-c(1:5)]
  
  assign("loci", loci, envir = .GlobalEnv)#Assigning loci to global environment so GCLr::combine_conflicts() can find it. 
  
  path <- paste0(QC_directory, "/Conflict Reports")
  
  qc_concordance_report_file <- list.files(path = path, pattern = "Concordance", full.names = TRUE)
  
  if(length(qc_concordance_report_file) == 0) {
    
    stop(paste0("No concordance files were found in: ", path, "\nPlease check to make sure that the directory path is correct and concordance files exist in the 'Conflict Reports' folder."))
    
  }
  
  for (silly in paste0(project_sillys,".gcl")) {
    
    # Create a new object name by removing ".gcl" and appending "_raw.gcl"
    new_silly_name <- paste0(stringr::str_remove(silly, ".gcl"), "_raw.gcl")
    
    # Assign the value of each silly.gcl to silly_raw.gcl
    assign(new_silly_name, get(silly), envir = .GlobalEnv)
    
  }
  
  conflicts <- GCLr::combine_conflicts(files = qc_concordance_report_file)
  
  success <- genotypes %>% 
    dplyr::select(plate_id = DNA_PLATE_ID, dplyr::all_of(loci)) %>% 
    tidyr::pivot_longer(-plate_id, names_to = "locus", values_to = "call") %>% 
    dplyr::group_by(plate_id) %>% 
    dplyr::summarize(success_rate = sum(!call %in% c("0", "0/0"))/length(call), .groups = "drop") %>% 
    dplyr::mutate(plate_id = as.character(plate_id))
  
  conf_succ <- conflicts %>% 
    dplyr::left_join(success, by = "plate_id") %>% 
    dplyr::rename(PlateID = plate_id)
  
  pivot1 <- conf_succ %>% 
    dplyr::filter(concordance_type %in% c("Homo-Het", "Het-Homo", "Homo-Homo", "Het-Het")) %>% 
    dplyr::select(SillySource, PlateID, success_rate, concordance_type) %>% 
    dplyr::group_by(SillySource, PlateID, success_rate, concordance_type) %>% 
    dplyr::summarize(conflict = length(concordance_type), .groups = "drop") %>% 
    dplyr::group_by(SillySource, PlateID, success_rate, concordance_type, conflict) %>% 
    dplyr::mutate(Total_conflicts = sum(conflict)) %>% 
    tidyr::pivot_wider(names_from = concordance_type, values_from = conflict) %>% 
    dplyr::select(-Total_conflicts, Total_conflicts)
  
  pivot2 <- conf_succ %>% 
    dplyr::filter(concordance_type %in% c("Homo-Het", "Het-Homo", "Homo-Homo", "Het-Het")) %>% 
    dplyr::select(locus, PlateID, success_rate, concordance_type) %>% 
    dplyr::group_by(locus, PlateID, success_rate, concordance_type) %>% 
    dplyr::summarize(conflict = length(concordance_type), .groups = "drop") %>% 
    dplyr::group_by(locus, PlateID, success_rate, concordance_type, conflict) %>% 
    dplyr::mutate(Total_conflicts = sum(conflict)) %>% 
    tidyr::pivot_wider(names_from = concordance_type, values_from = conflict) %>% 
    dplyr::select(-Total_conflicts, Total_conflicts)
  
  #Need to append the two pivot tables to the genotypes report xlsx file
  xlfile <- paste0(project, "_genotypes_report.xlsx")
  
  wb <- xlsx::loadWorkbook(xlfile)
  
  bouncexlfile <- paste0(project, "QC_bounce.xlsx")
  
  xlsx::saveWorkbook(wb, file = bouncexlfile)
  
  xlsx::write.xlsx(x = pivot1 %>% as.data.frame(), file = bouncexlfile, append = TRUE, sheetName = "Conflicts by individual", showNA = FALSE, row.names = FALSE)
  xlsx::write.xlsx(x = pivot2 %>% as.data.frame(), file = bouncexlfile, append = TRUE, sheetName = "Conflicts by locus", showNA = FALSE, row.names = FALSE)
  
}