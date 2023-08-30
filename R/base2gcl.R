#' Create `.gcl` Objects from a \emph{rubias} Baseline Object
#' 
#' This function creates a `.gcl` object for each collection in a \emph{rubias} baseline object.
#' 
#' @param base a rubias baseline (aka reference) object.
#' 
#' @return `.gcl` objects are assigned to your global environment and a vector of silly codes (aka `sillyvec`) is returned (see details)
#' 
#' @details
#' This function is a quick way of creating `.gcl` objects from a \emph{rubias} baseline object. The objects will contain the same variables as `.gcl` objects produced by [GCLr::loki2r] and [GCLr::loki2r_gaps] for compatibility with the [GCLr] package; however, most of the attribute variables will contain NAs.  
#' 
#' @examples
#' 
#' GCLr::base2gcl(base = GCLr::ex_baseline)
#'
#' @export
base2gcl <- function(base){
  
  #Check if supplied object
  vars.check <- base %>% 
    names()
  
  if(!sum(vars.check[1:4] %in% c("sample_type","collection","indiv")) == 3){
    
    stop("The baseline object must contain the following variables: sample_type, collection, indiv")
    
  }
  
  if(!base$sample_type %>% unique() == "reference"){
    
    stop("The supplied object does not appear to be a baseline (aka reference) object.")
    
  }
  
  #These are the 19 attribute variables that all .gcl objects must have.
  
  attr.vars <- c("FK_FISH_ID", "COLLECTION_ID", "SILLY_CODE", "PLATE_ID", "PK_TISSUE_TYPE", 
                 "CAPTURE_LOCATION", "CAPTURE_DATE", "END_CAPTURE_DATE", "MESH_SIZE", 
                 "MESH_SIZE_COMMENT", "LATITUDE", "LONGITUDE", "AGENCY", "VIAL_BARCODE", 
                 "DNA_TRAY_CODE", "DNA_TRAY_WELL_CODE", "DNA_TRAY_WELL_POS", "CONTAINER_ARRAY_TYPE_ID", 
                 "SillySource") 
  
  #Locus variables
  loc.vars <- base %>% 
    dplyr::select(-sample_type, -repunit, -collection, -indiv) %>% 
    names()
  
  all.df <- base %>% 
    dplyr::mutate(SillySource = indiv,
                  SILLY_CODE = collection) %>% 
    tidyr::separate(indiv, into = c(NA, "FK_FISH_ID")) %>% 
    dplyr::group_by(collection) %>% 
    dplyr::mutate(FK_FISH_ID = seq(length(FK_FISH_ID)),
                  COLLECTION_ID = NA_real_,
                  PLATE_ID = NA_character_,
                  PK_TISSUE_TYPE = NA_character_,
                  CAPTURE_LOCATION = NA_character_,
                  CAPTURE_DATE = NA,
                  END_CAPTURE_DATE = NA,
                  MESH_SIZE = NA_character_,
                  MESH_SIZE_COMMENT = NA_character_,
                  LATITUDE = NA_real_,
                  LONGITUDE = NA_real_,
                  AGENCY = NA_character_,
                  VIAL_BARCODE = NA_character_,
                  DNA_TRAY_CODE = NA_character_,
                  DNA_TRAY_WELL_CODE = NA_real_,
                  DNA_TRAY_WELL_POS = NA_character_,
                  CONTAINER_ARRAY_TYPE_ID = NA_real_) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(dplyr::all_of(c(attr.vars, loc.vars)))
  
  sillys <- all.df$SILLY_CODE %>% unique()
  
  for(silly in sillys){
    
    assign(x = paste0(silly, ".gcl"), value = all.df %>% dplyr::filter(SILLY_CODE == silly), pos = -1, envir = .GlobalEnv)
    
  }
  
  print(paste0("The following '.gcl' objects have been assigned to your global environment: ", paste0(sillys, ".gcl", collapse = ", ")))
  
  return(sillys)
  
}