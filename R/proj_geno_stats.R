#' Calculate Project Genotype Statistics and Write to Excel Workbook
#' 
#' This function takes the genotypes report object from `GCLr::get_geno()`, calculates genotypes statistics, and write them to a formatted Excel workbook.
#' 
#' @param report a report object created by `GCLr::get_geno()`
#' @param file the name of the excel workbook.
#' 
#' @details This function is called on by `GCLr::get_geno()` when a `project_name` is supplied and `proj.stats = TRUE`.  However, this function can also be run separately if you have an genotypes object created by `GCLr::get_geno`.
#' 
#' @seealso [GCLr::get_geno()]
#' 
#' @returns This function writes an Excel workbook with three sheets: 
#' \itemize{
#'   \item \code{proj_Genotypes_Table}: a formatted genotypes table
#'   \item \code{Calculations}: number of fish with success rates < 80%, total number of project fish, and Proportion of fish with success rates greater than 80%
#'   \item \code{Avg success by plate}: a table with PlateID and Average Success Rate
#'   }
#' 
#' @examples
#' \dontrun{
#' Geno_proj <- GCLr::get_geno(project_name = "K205", file = NULL, username = "awbarclay", password = "password", open.file = TRUE) 
#' 
#' GCLr::proj_geno_stats(report = Geno_proj, file = path.expand("~/TestGenotypesReport_proj.xlsx"))
#' }   
#'   
#' @export
proj_geno_stats <- function(report, file){
  
 
   #Checking to make sure file has the .xlsx extension
  
  if(grep(pattern = "*.xlsx", x = file)  %>% length() == 0){
    
    stop("This function creates an Excel file. Please supply an file path with a .xlsx extension")
    
  }
  
  if(sum(grepl("DNA_PLATE_ID", names(report))) == 0){
    
    stop("The genotypes report supplied does not contain DNA_PLATE_ID. 
         Was the report produced by GenotypeReport.GCL using the project_name argument?")
    
  }
  
  # Excel Sheet 1 - Add count and success rate variables
  
  project <- report$LAB_PROJECT_NAME %>% unique()
  
  sheet1_name = paste0(paste0(project, collapse = "_"), "_Genotypes_Table")
  
  loc_vars <- names(report)[-c(1:5)]# This is for selecting columns in the correct order
  
  report_out <- report %>%  
    tidyr::pivot_longer(dplyr::all_of(loc_vars), names_to = "locus", values_to = "geno") %>% 
    dplyr::group_by(LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, DNA_PLATE_ID) %>% 
    dplyr::mutate(`Count of 0/0` = stringr::str_count(geno,"0/0") %>% sum()) %>% 
    tidyr::pivot_wider(names_from = locus, values_from = geno) %>% 
    dplyr::mutate(`Success Rate` = 1-(`Count of 0/0`/length(loc_vars)), SILLY_CODE_FK_FISH_ID = paste(SILLY_CODE, FK_FISH_ID, sep = "_"), Blank = ".") %>% 
    dplyr::select(LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, FK_FISH_ID, SILLY_CODE_FK_FISH_ID, PlateID = DNA_PLATE_ID, `Count of 0/0`, `Success Rate`, Blank, dplyr::all_of(loc_vars))
  
  # Excel Sheet 2 - calculations
  sheet2_name = "Calculations"
  
  below80 <- sum(report_out$`Success Rate` < 0.8)
  
  total_fish <- length(report_out$`Success Rate`)
  
  overall_success <- 1-(below80/total_fish) 
  
  calcs_out <- tibble::tibble(stat_name = c("less than 80%", "total fish", "Proportion of Fish > 80%"), stat = c(below80, total_fish, overall_success))
  
  # Excel Sheet 3 - Average Success Rate by plate ID
  
  sheet3_name <- "Avg success by plate"
  
  avg_success <- report_out %>% 
    dplyr::group_by(PlateID) %>% 
    dplyr::summarise(`Average Success Rate` = mean(`Success Rate`))
  
  # Create workbook
  
  wb <- openxlsx::createWorkbook() 
  
  # Add sheet 1
  
  openxlsx::addWorksheet(wb, sheetName = sheet1_name ) 
  
  openxlsx::writeData(wb, sheet = sheet1_name, x = report_out, startCol = 1, startRow = 1, colNames = TRUE)
  
  head_style <-  openxlsx::createStyle(textRotation = 90, halign = "center")
  
  hallign <-  openxlsx::createStyle(halign = "center")
  
  blankcol <-  openxlsx::createStyle(textRotation = 90, fontColour = "#D9D9D9", fgFill = "#D9D9D9")
  
  openxlsx::addStyle(wb, sheet = sheet1_name, style = head_style, rows = 1, cols = 1:ncol(report_out))
  
  openxlsx::addStyle(wb, sheet = sheet1_name, style = hallign, rows = 2:(nrow(report_out)+1), cols = 10:ncol(report_out), gridExpand = TRUE) 
  
  openxlsx::addStyle(wb, sheet = sheet1_name, style = blankcol, rows = 1:(nrow(report_out)+1), cols = 9, gridExpand = TRUE) 
  
  # Add sheet 2
  
  openxlsx::addWorksheet(wb, sheetName = sheet2_name ) 
  
  openxlsx::writeData(wb, sheet = sheet2_name, x = calcs_out, startCol = 1, startRow = 1, colNames = FALSE)
  
  openxlsx::conditionalFormatting(wb, sheet = sheet1_name, cols = 8, rows = 1:total_fish, rule = "< .8", style =  openxlsx::createStyle(fontColour = "#806000", bgFill = "#FFE699"))# conditional formatting for Success rate below 80%
  
  openxlsx::conditionalFormatting(wb, sheet = sheet1_name, cols = 1:length(loc_vars)+9, rows = 2:(total_fish+1), type = 'contains', rule = "0/0") # conditional formatting for 0/0s in red - catches markers that get 0's
  
  openxlsx::conditionalFormatting(wb, sheet = sheet1_name, cols = 1:length(loc_vars)+9, rows = 2:(total_fish+1), type = 'contains', rule = "0") # conditional formatting for 0s in red - this catches mt markers that are 0
  
  # Add sheet 3
  
  openxlsx::addWorksheet(wb, sheetName = sheet3_name)
  
  openxlsx::writeData(wb, sheet = sheet3_name, x = avg_success , startCol = 1, startRow = 1, colNames = TRUE)
  
  # Set column widths for sheet1
  
  width_vec <- apply(as.data.frame(report_out), 2, function(x) max((nchar(as.character(x), allowNA = TRUE, keepNA = NA) ), na.rm = TRUE))#Get max width of characters is each column
  
  openxlsx::setColWidths(wb, sheet = sheet1_name, cols = 1:ncol(report_out), widths = width_vec)
  
  # Save workbook
  
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  
} 