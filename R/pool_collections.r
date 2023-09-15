#' Combine "*.gcl" Objects
#'
#' This function combines "*.gcl" objects into a new one called "newname.gcl".
#'
#' @param collections A character vector of silly codes without the ".gcl" extension.
#' 
#' @param loci A character vector of locus names.
#' 
#' @param IDs A named list of FK_FISH_ID
#'             These will be used to subset each collection before pooling. If no IDs are supplied, all individuals from each collection are used.
#'             
#' @param newname The name of the new "*.gcl" created. Do not provide the ".gcl" extension. If no name is supplied, then the newname defaults to
#'               the collection names collapsed with a period between each name (e.g., "KQUART06.KQUART08.KQUART09").
#'
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)  
#'
#' @details
#' This function is useful for pooling collections from the same location into a single population ("pooled collection") and for producing "pooled mixture" objects for mixed stock analysis.
#'
#' @return assigns a new "pooled collection" to your global environment
#' 
#' @examples
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline, unpool = TRUE)
#' 
#' loci <- GCLr::ex_LocusControl$locusnames[-c(10, 12, 13, 32, 33)]
#' 
#' GCLr::pool_collections(collections = sillyvec[1:2], loci = loci, IDs = NULL, newname =  paste(sillyvec[1:2], collapse = "."), LocusCtl = GCLr::ex_LocusControl)
#' 
#' @export
pool_collections <- function(collections, loci = LocusControl$locusnames, IDs = NULL, newname = paste(collections, collapse = "."), LocusCtl = LocusControl){
  
  if(nchar(newname) > 200){
    
    newname <- substr(newname, start = 1, stop = 200)
    
  }
  
  if(!all(loci %in% LocusCtl$locusnames)){
    
    stop(paste0("The following `loci` were not found in `LocusControl`:\n", paste(setdiff(loci, LocusCtl$locusnames), collapse = "\n")))
    
  }
  
  ncollections <- length(collections)
  
  if(is.null(IDs)){
    
    IDs <- sapply(collections, function(collection){
      
      get(paste0(collection, ".gcl"), pos = 1)$FK_FISH_ID
      
    }, simplify = FALSE) 
    
  }
  
  if(!is.list(IDs)){
    
    stop("'IDs' must be a list")
    
  }
  
  if(ncollections != length(IDs)){
    
    stop("'IDs' must be same length as 'collections'")
    
  }
  
  IDs <- purrr::set_names(IDs, collections)  # Making sure IDs has names
  
  SubsetLoci <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% as.vector()  # These are the locus score headers for subsetting by loci.
  
  output <- lapply(collections, function(collection){
    
    my.gcl <- get(paste0(collection, ".gcl"), pos = 1)
    
    attr <-  c("FK_FISH_ID", "COLLECTION_ID", "SILLY_CODE", "PLATE_ID", "PK_TISSUE_TYPE", "CAPTURE_LOCATION", "CAPTURE_DATE", "END_CAPTURE_DATE", "MESH_SIZE", "MESH_SIZE_COMMENT", "LATITUDE", "LONGITUDE", "AGENCY", "VIAL_BARCODE", "DNA_TRAY_CODE", "DNA_TRAY_WELL_CODE", "DNA_TRAY_WELL_POS", "CONTAINER_ARRAY_TYPE_ID", "SillySource")  # required attribute names
    
    my.IDs <- IDs[[collection]]
    
    my.gcl %>% 
      dplyr::filter(FK_FISH_ID %in% my.IDs) %>% 
      dplyr::select(tidyselect::all_of(attr), tidyselect::all_of(SubsetLoci))
    
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(FK_FISH_ID = seq(length(unlist(IDs))), SILLY_CODE = newname)
  
  assign(paste0(newname, ".gcl"), output, pos = 1, envir = .GlobalEnv) 
  
  message(paste0(newname, ".gcl has been assigned to your global environment."))
  
}