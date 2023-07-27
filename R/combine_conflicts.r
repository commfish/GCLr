#' @title Combine QC Conflict Reports and Add PLATE_ID Column
#' @description This function combines QC conflict reports and adds a column of PLATE_IDs associated with each individual conflict.  The function requires "*.gcl" objects and the "ProjectSillys" vector created by ReadProjectloki2r.
#'
#' @param files A vector of CSV files that you want to combine, including the file extension.
#'
#' 
#' @family QA scripts
#' 
#' @examples
#'  \dontrun{
#'  
#' combine_conflicts(files = c("CM31 qc 01 Import Results.csv", "CM31 qc 02 Import Results.csv", "CM31 qc Plate3 Imports.csv"))
#' }
#' 
#' @export
combine_conflicts <- function(files) {
  
  if (sys.call()[[1]] == quote(CombineConflictsWithPlateID.GCL.r)) {
    warning("The function 'CombineConflictsWithPlateID.GCL.r' is deprecated. Please use 'combine_conflicts' instead.")
  }
  # Read in concordance files
  suppressMessages(
    concordance <- files %>% 
      purrr::map(readr::read_csv) %>% 
      dplyr::bind_rows()
  )
  
  # Rename headers if new format
  if ( "Sample Number" %in% names(concordance)) {
    concordance <- concordance %>% 
      dplyr::rename("Fish ID" = "Sample Number", 
                    "File Allele 1" = "File: Allele 1",
                    "File Allele 2" = "File: Allele 2",
                    "DB Allele 1" = "Database: Allele 1",
                    "DB Allele 2" = "Database: Allele 2"
      )
  }
  
  # Re-format concordance tibble
  concordance <- concordance %>% 
    tidyr::unite(silly_source, c(`Silly Code`, `Fish ID`), sep = "_", remove = FALSE) %>% 
    dplyr::mutate(`Silly Code` = factor(x = `Silly Code`, levels = ProjectSillys))
  
  # Old conflict report has "0" for mitochondrial conflicts, new has " " for mitochondrial conflicts, we will refer to them as "Homo-Homo".
  level_key <- list(`DB Zero` = "DB Zero", `File Zero` = "File Zero", `Het-Het` = "Het-Het", `Het-Homo` = "Het-Homo", `Homo-Het` = "Homo-Het", `Homo-Homo` = "Homo-Homo", `0` = "Homo-Homo", ` ` = "Homo-Homo")
  types <- c("DB Zero", "File Zero", "Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo")  # order with levels
  
  # Pool all collections in to one master silly
  pool_collections(collections = ProjectSillys, loci = loci, newname = "master")
  
  # Tibble of reduced attributes table
  master.tbl <- master.gcl$attributes %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(SillySource, PLATE_ID)
  
  # Join plate id
  concordance_join <- concordance %>% 
    dplyr::left_join(master.tbl, by = c("silly_source" = "SillySource")) %>% 
    dplyr::rename(silly = `Silly Code`, 
                  fish_id = `Fish ID`,
                  locus = Locus,
                  file_allele_1 = `File Allele 1`,
                  file_allele_2 = `File Allele 2`,
                  db_allele_1 = `DB Allele 1`,
                  db_allele_2 = `DB Allele 2`,
                  concordance = Concordance,
                  concordance_type = `Concordance Type`,
                  plate_id = PLATE_ID) %>% 
    dplyr::select(silly, fish_id, silly_source, locus, file_allele_1, file_allele_2, db_allele_1, db_allele_2, concordance, concordance_type, plate_id) %>% 
    dplyr::filter(concordance == "Conflict") %>% 
    dplyr::mutate(concordance_type = dplyr::recode_factor(concordance_type, !!!level_key)) %>%  # recode to deal with mitochondrial conflicts
    dplyr::mutate(concordance_type = factor(x = concordance_type, levels = types)) %>%  # new levels
    dplyr::mutate(locus = factor(x = locus, levels = loci)) %>%  # new levels
    dplyr::mutate(silly = factor(x = silly, levels = ProjectSillys)) %>% 
    dplyr::mutate(plate_id = factor(x = plate_id, levels = sort(unique(plate_id))))
  
  # Write copy
  readr::write_csv(concordance_join, path = "Conflict Reports/CombinedConflictsWithPlateID.csv")
  
  # Assign to global environment
  assign(x = "combined_conflicts", value = concordance_join, pos = 1)
  
}