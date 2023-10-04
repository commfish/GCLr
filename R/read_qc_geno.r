#' @title Read QC Data
#'
#' @description This function reads QC data from Biomark CSV files. It's called on by [GCLr::qc()].
#'
#' @param qc_csv_filepaths character vector specifying the file paths of the QC CSV files.
#' @param skip number of lines to skip while reading the CSV files when type = "Biomark". (default = 15).
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)  
#' @param type the type of project ("Biomark", "uSat", or "GT-seq"), (default = "Biomark") 
#'
#' @returns Returns a few silly objects to the global environment 
#'   - `qc.gcl` tibble objects
#'   - `qc_sillys`; a character vector of qc sillys
#'
#' @export
read_qc_geno <- function(qc_csv_filepaths, skip = 15, LocusCtl = LocusControl, type = c("Biomark", "uSat", "GT-seq")[1]) {
  
  sillyvec <- as.vector(sapply(objects(pattern = "*\\.gcl", pos = 1), function(gclname){strsplit(gclname, split = "\\.gcl")[[1]][1]}))
  
  if(type == "Biomark"){
    
    Genotypesqc <- lapply(qc_csv_filepaths, function(pth){
      
      read.table(pth, header = TRUE, sep = ",", colClasses = "character", stringsAsFactors = FALSE, skip = skip)[ , c("Name", "Converted", "Assay", "ID", "Allele.X", "Allele.Y")]
      
    }) %>% dplyr::bind_rows() %>% 
      dplyr::filter(Name != "NTC") %>%
      dplyr::mutate(Converted = dplyr::case_when(Converted == "No Call"| Converted == "Invalid" | Converted == "0" ~ NA_character_,
                                                 TRUE~Converted)) %>% 
      tidyr::separate(Name, into = c("SILLY_CODE", "FK_FISH_ID"), sep = "_", remove = TRUE) %>% 
      dplyr::mutate(FK_FISH_ID = as.numeric(FK_FISH_ID)) %>% 
      dplyr::mutate(SILLY_CODE = paste0(SILLY_CODE, "qc")) %>%  
      dplyr::mutate(SillySource = paste0(SILLY_CODE, "_", FK_FISH_ID)) %>% 
      tidyr::separate(Converted, into = c("allele1", "allele2"), sep = ":") %>% 
      dplyr::select(SillySource, SILLY_CODE, FK_FISH_ID, allele1, allele2, locus = Assay) 
    
  }
  
  if(type == "uSat"){
    
    Genotypesqc <- lapply(qc_csv_filepaths, function(pth){
      
      suppressMessages(readr::read_csv(file = pth, show_col_types = FALSE)[, 1:4])
      
    }) %>% dplyr::bind_rows() %>% 
      dplyr::rename(locus = Marker) %>% 
      dplyr::mutate(allele1 = as.character(`Allele 1`), allele2 = as.character(`Allele 2`)) %>% 
      tidyr::separate(`Sample Name`, into = c("SILLY_CODE", "FK_FISH_ID"), sep = "_", remove = FALSE) %>% 
      dplyr::mutate(SILLY_CODE = paste0(SILLY_CODE, "qc")) %>% 
      dplyr::mutate(SillySource = paste0(SILLY_CODE, "_", FK_FISH_ID), FK_FISH_ID = as.numeric(FK_FISH_ID)) %>% 
      dplyr::select(SillySource, SILLY_CODE, FK_FISH_ID, allele1, allele2, locus) 
  }
  
  if(type == "GT-seq"){
    
    Genotypesqc <- lapply(qc_csv_filepaths, function(pth){
      
      suppressMessages(readr::read_csv(file = pth, show_col_types = FALSE, na = c("", "NA", "0", "0/0")))
      
    }) %>% dplyr::bind_rows() %>% 
      dplyr::rename(FK_FISH_ID = SAMPLE_NUM, locus = LOCUS) %>% 
      dplyr::mutate(SILLY_CODE = paste0(SILLY_CODE, "qc")) %>% 
      tidyr::unite(SillySource, c(SILLY_CODE, FK_FISH_ID), sep = "_", remove = FALSE) %>% 
      tidyr::separate(GENOTYPE, into = c("allele1", "allele2"), sep = "/") %>% 
      dplyr::select(SILLY_CODE, FK_FISH_ID, SillySource, locus, allele1, allele2) 
  }
  
  # Add extra columns as NAs for compatibility with all other functions
  attnames <- c("FK_FISH_ID", "COLLECTION_ID", "SILLY_CODE", "PLATE_ID", "PK_TISSUE_TYPE", 
                "CAPTURE_LOCATION", "CAPTURE_DATE", "END_CAPTURE_DATE", "MESH_SIZE", 
                "MESH_SIZE_COMMENT", "LATITUDE", "LONGITUDE", "AGENCY", "VIAL_BARCODE", 
                "DNA_TRAY_CODE", "DNA_TRAY_WELL_CODE", "DNA_TRAY_WELL_POS", "CONTAINER_ARRAY_TYPE_ID", 
                "SillySource") # list of 19 allowed attributes
  
  sillyvecqc0 <- unique(sub(x = Genotypesqc$SILLY_CODE, pattern = "*qc", replacement = ""))
  
  sillyvecqc <- paste0(sillyvecqc0[!is.na(match(sillyvecqc0, sillyvec))], "qc")
  
  for(silly in sillyvecqc){
    
    my.dat <- Genotypesqc %>% 
      dplyr::filter(SILLY_CODE == silly) %>% 
      tidyr::pivot_wider(names_from = locus, values_from = c("allele1", "allele2"), names_glue = "{locus}_{.value}", values_fill = NULL) %>% 
      dplyr::rename_with(~ gsub("_allele1", "", .), dplyr::ends_with("_allele1")) %>% #fix names locus
      dplyr::rename_with(~ gsub("_allele2", ".1", .), dplyr::ends_with("_allele2")) #fix names locus.1
    
    # Check which attnames exist in the data_master tibble, identify missing
    missing_columns <- setdiff(attnames, names(my.dat))
    
    # Add the missing attnames, as NA values, and arrange (19 attributes first)
    my.dat <- my.dat %>%
      tibble::add_column(.data = ., !!!setNames(rep(NA, length(missing_columns)), missing_columns)) %>% 
      dplyr::select(dplyr::all_of(attnames), sort(names(.), na.last = TRUE))
    
    assign(paste(silly, ".gcl", sep = ""), my.dat, pos = 1)
    
    message(paste0(silly, ".gcl created ", match(silly, sillyvecqc), " of ", length(sillyvecqc), " completed."))
    
    }
  
  assign(x = "qc_sillys", value = sillyvecqc, pos = 1)
  
} 