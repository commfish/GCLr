old2new_gcl <- function(sillyvec, save_old = FALSE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function converts old style "*.gcl" objects from nested arrays to tibbles of scores and attributes. 
  # The old style objects can be saved as *.gcl_old before they are overwritten.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - A character vector of silly codes with associated *.gcl objects that need to be converted.
  #   overwrite - logical; whether you want the old objects overwritten without saving 
  #              (TRUE) or assign the new objects to *.gcl_old before overwriting (FALSE) 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function assigns a converted *.gcl object to the current workspace.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  # 
  # old2new_gcl(sillyvec = sillyvec157, save_old = TRUE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  if(!all(sillyvec %in% stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))) {  # Do all sillys exist in the environment?
    
    missing_sillys <- setdiff(sillyvec, stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))
    
    stop(paste0("The following sillys are not in your environment:\n", paste0(missing_sillys, collapse = ", ")))
    
  }
  
  sillyvec0 <- sillyvec[sapply(sillyvec, function(silly){  # Excluding objects that are already in tidy format.
    
    !tibble::is_tibble(get(paste0(silly, ".gcl")))
    
  })]
  
  attr_sel <- c("FK_FISH_ID", "COLLECTION_ID", "SILLY_CODE", "PLATE_ID", "PK_TISSUE_TYPE", 
    "CAPTURE_LOCATION", "CAPTURE_DATE", "END_CAPTURE_DATE", "MESH_SIZE", 
    "MESH_SIZE_COMMENT", "LATITUDE", "LONGITUDE", "AGENCY", "VIAL_BARCODE", 
    "DNA_TRAY_CODE", "DNA_TRAY_WELL_CODE", "DNA_TRAY_WELL_POS", "CONTAINER_ARRAY_TYPE_ID", 
    "SillySource") #These are the 19 attribute fields that are allowed.
  
  sapply(sillyvec0, function(silly){
    
    my.gcl <- get(paste0(silly, ".gcl"))
    
    if(save_old){assign(paste0(silly, ".gcl_old"), value = my.gcl, pos = -1, envir = .GlobalEnv)}  # Saving old gcl
    
    my_attr <- my.gcl$attributes %>% 
      names()
    
    if(sum(!attr_sel %in% my_attr) > 0){ # Rename some columns if they are different and add missing columns.
      
      if("plateID" %in% my_attr){
        
        my.gcl$attributes <- my.gcl$attributes %>% 
          dplyr::rename(PLATE_ID = plateID)
        
      }
      
      if("FK_COLLECTION_ID" %in% my_attr){
        
        my.gcl$attributes <- my.gcl$attributes %>% 
          dplyr::rename(COLLECTION_ID = FK_COLLECTION_ID)
        
      }
      
      if(!"CONTAINER_ARRAY_TYPE_ID" %in% my_attr){
       
        my.gcl$attributes$CONTAINER_ARRAY_TYPE_ID <- NA_real_
        
      }
      
      if(!"LATITUDE" %in% my_attr){
        
        my.gcl$attributes$LATITUDE <- NA_real_
        
      }
      
      if(!"LONGITUDE" %in% my_attr){
        
        my.gcl$attributes$LONGITUDE <- NA_real_
        
      }
      
      if(!"AGENCY" %in% my_attr){
        
        my.gcl$attributes$AGENCY <- NA_character_
        
      }
      
      if(!"DNA_TRAY_CODE" %in% my_attr){
        
        my.gcl$attributes$DNA_TRAY_CODE <- NA_character_
        
      }
      
      if(!"DNA_TRAY_WELL_CODE" %in% my_attr){
        
        my.gcl$attributes$DNA_TRAY_WELL_CODE <- NA_real_
        
      }
      
      if(!"DNA_TRAY_WELL_POS" %in% my_attr){
        
        my.gcl$attributes$DNA_TRAY_WELL_POS <- NA_character_
        
      }
      
      if(!"MESH_SIZE_COMMENT" %in% my_attr){
        
        my.gcl$attributes$MESH_SIZE_COMMENT <- NA_character_
        
      }
      
    }
    
    attr <- my.gcl$attributes %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(FK_FISH_ID = as.double(FK_FISH_ID),  # FK_FISH_ID must be numeric
                    COLLECTION_ID = as.double(COLLECTION_ID),
                    SILLY_CODE = silly,  # Some older attributes did not include SILLY_CODE
                    PLATE_ID = as.character(PLATE_ID),
                    PK_TISSUE_TYPE = as.character(PK_TISSUE_TYPE),
                    CAPTURE_LOCATION = as.character(CAPTURE_LOCATION),
                    CAPTURE_DATE = lubridate::as_date(CAPTURE_DATE),  # make date
                    END_CAPTURE_DATE = lubridate::as_date(END_CAPTURE_DATE),
                    MESH_SIZE = as.character(MESH_SIZE),
                    MESH_SIZE_COMMENT = as.character(MESH_SIZE_COMMENT),
                    LATITUDE = as.double(LATITUDE),
                    LONGITUDE = as.double(LONGITUDE),
                    AGENCY = as.character(AGENCY),
                    VIAL_BARCODE = as.character(VIAL_BARCODE),
                    DNA_TRAY_CODE = as.character(DNA_TRAY_CODE),
                    DNA_TRAY_WELL_CODE = as.double(DNA_TRAY_WELL_CODE),
                    DNA_TRAY_WELL_POS = as.character(DNA_TRAY_WELL_POS),
                    CONTAINER_ARRAY_TYPE_ID = as.double(CONTAINER_ARRAY_TYPE_ID),
                    SillySource = as.character(SillySource)) %>% 
      dplyr::select(all_of(attr_sel)) # Select the 19 attribute columns
    
    scores <- lapply(seq(dim(my.gcl$scores)[[3]]), function(dim){
      
      s = my.gcl$scores[ , , dim, drop = FALSE] %>%
        tibble::as_tibble()
      
      colnames(s) <- dimnames(my.gcl$scores)[[2]] 
      
      if(dim > 1){
        
        colnames(s) <- paste(colnames(s), dim - 1, sep = ".")
        
      }
      
      s
      
    }) %>% 
      dplyr::bind_cols() %>% 
      dplyr::mutate(dplyr::across(dplyr::where(is.character), ~dplyr::na_if(., "0")))
    
    tidy.gcl <- scores %>%
      dplyr::bind_cols(attr) %>% 
      dplyr::select(colnames(attr), sort(colnames(scores)))
    
    assign(paste0(silly, ".gcl"), value = tidy.gcl, pos = -1, envir = .GlobalEnv)
    
  })  # End silly
  
  message(paste0("The following *.gcl objects have been converted to tibbles:\n", paste0(sillyvec, collapse = ", ")))
  
}