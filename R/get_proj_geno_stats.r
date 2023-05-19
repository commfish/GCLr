get_proj_geno_stats <- function(report = NULL, in.file = NULL, out.file = "GenotypeReportqc.xlsx"){  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function takes an unmodified genotypes report object or file produced by get_geno and
  #  calculates statistics using by GCL lab staff for checking project genotypes, then writes out a formatted .xlsx file of the stats.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #  report - a project genotypes report object created by get_geno (i.e., genotypes pulled using project name)
  #   
  #  in.file - the file path to an unmodified genotypes report .csv file created by get_geno
  #
  #  out.file - the file path for the .xlsx output file with file name and extension.
  #             If not supplied, the file named "GenotypeReportqc.xlsx" will be written to your current working directory
  #             Check your current directory with getwd()
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #  xlsx file with three sheets
  # 
  #  Sheet1 - the report produced by get_geno when pulling by project_name with "Count of 0/0" and "Success Rate" columns added and conditional formatting.
  #  
  #  Sheet2 - a.k.a "Calculations"
  # 
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # .password <- "************"
  # username <- "awbarclay"
  # project_name <- "P061"
  #
  # GenotypeReport <- get_geno(sillyvec = NULL, loci = NULL, path = "GenotypeReport_test.csv", username = username, password = .password, project_name = project_name, open.file = FALSE)
  #
  # get_proj_geno_stats(report = GenotypeReport, in.file = NULL, out.file = "FormattedGeneotypesReport_test.xlsx")
  #
  # get_proj_geno_stats(report = NULL, in.file = "GenotypeReport.csv", out.file = "FormattedGeneotypesReport_test.xlsx")
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Checking that arguments are supplied.
  if(is.null(report) & is.null(in.file)){ 
    
    stop("The user must supply a report object or the file path to a genotypes report .csv file.")
    
  }
  
  # Checking that only one input is supplied
  if(!is.null(report) & !is.null(in.file)){ 
    
    stop("Both report and in.file arguments were supplied.  The user must supply one or the other, not both.")
    
  }
  
  # Requiring the tidyverse and openxlsx packages. Packages will be installed if the user doesn't have them installed
  
  #Checking to make sure the out.file has the .xlsx extension
  
  if(grep(pattern = "*.xlsx", x = out.file)  %>% length() == 0){
    
    stop("Please supply an out.file with a .xlsx extension")
    
  }
    
  # If in.file is supplied read in report
  if (!is.null(in.file) & is.null(report)) {
    
    report <- readr::read_csv(in.file, col_types = cols(
      LAB_PROJECT_NAME = col_character(),
      FK_COLLECTION_ID = col_double(),
      SILLY_CODE = col_character(),
      FK_FISH_ID = col_double(),
      .default = "c")
      )
  }
  

  if(sum(grepl("DNA_PLATE_ID", names(report))) == 0){
    
    stop("The genotypes report supplied does not contain DNA_PLATE_ID. 
         Was the report produced by get_geno using the project_name argument?")
    
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
    dplyr::select(LAB_PROJECT_NAME, FK_COLLECTION_ID, SILLY_CODE, SILLY_CODE_FK_FISH_ID, FK_FISH_ID, PlateID = DNA_PLATE_ID, `Count of 0/0`, `Success Rate`, Blank, dplyr::all_of(loc_vars))
  
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
  
  openxlsx::conditionalFormatting(wb, sheet = sheet1_name, cols = 8, rows = 1:total_fish, rule = "< .8", style =  openxlsx::createStyle(fontColour = "#806000", bgFill = "#FFE699"))
  
  openxlsx::conditionalFormatting(wb, sheet = sheet1_name, cols = 1:length(loc_vars)+9, rows = 2:(total_fish+1), type = 'contains', rule = "0/0") # conditional formatting for 0/0s in red - catches markers that get 0's
  
  openxlsx::conditionalFormatting(wb, sheet = sheet1_name, cols = 1:length(loc_vars)+9, rows = 2:(total_fish+1), type = 'contains', rule = "0") # conditional formatting for 0s in red - this catches mt markers that are 0
  
  # Add sheet 3
  
  openxlsx::addWorksheet(wb, sheetName = sheet3_name)
  
  openxlsx::writeData(wb, sheet = sheet3_name, x = avg_success , startCol = 1, startRow = 1, colNames = TRUE)
  
  # Set column widths for sheet1
  
  width_vec <- apply(as.data.frame(report_out), 2, function(x) max((nchar(as.character(x), allowNA = TRUE, keepNA = NA) ), na.rm = TRUE))#Get max width of characters is each column
  
  openxlsx::setColWidths(wb, sheet = sheet1_name, cols = 1:ncol(report_out), widths = width_vec)
  
  # Save workbook
  
  openxlsx::saveWorkbook(wb, out.file, overwrite = TRUE)
  
}
  