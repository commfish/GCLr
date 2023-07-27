#' Remove IDs from the ".gcl" object associated with "silly"
#'
#' This function removes specific IDs from a Silly GCL object based on the provided 'IDs' vector.
#'
#' @param silly Character string specifying the name of the Silly GCL object. The function will attempt to find the object with the name 'silly.gcl'.
#' @param IDs Numeric or character vector containing the IDs to be removed from the Silly GCL object. If 'IDs' is a character vector, it will be converted to numeric.
#'
#' @return A tibble with information about the removal process.
#'   \item{SILLY_CODE}{Character string, the name of the Silly GCL object.}
#'   \item{IDs_Removed}{List, a list containing the IDs that were removed from the Silly GCL object. The list is named with the value of 'silly'.}
#'   \item{is_empty}{Logical, indicating whether the Silly GCL object is now empty after the removal.}
#'
#' @details
#' The 'remove_ids' function removes specific IDs from the .gcl object associated with 'silly'. The function checks if the 'silly' object is empty and raises a stop error if it is empty, as the removal process would not be possible.
#'
#' If 'IDs' contains character values, they will be converted to numeric. The function then checks if all the 'IDs' provided are present in the Silly GCL object. If any of the IDs are not found, it raises a stop error.
#'
#' After removing the specified IDs from the Silly GCL object, the function gives a message indicating the number of IDs removed and their respective 'silly' object name.
#'
#' Additionally, if the Silly GCL object becomes empty after the removal, a warning is issued to notify the user.
#'
#' @examples
#' \dontrun{
#' password <- "************"
#' 
#' create_locuscontrol(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
#' sillyvec <- c("SMCDO03", "SNEVA13")
#' loki2r(sillyvec = sillyvec, username = "awbarclay", password = password)
#' remove_ind_miss_loci(sillyvec = sillyvec)
#' 
#' remove_ids(silly = "SMCDO03", IDs = 1:10)
#' }
#'
#' @seealso
#' \code{\link{add_ids}}: Function to add IDs to a Silly GCL object.
#'
#' @export
  
remove_ids <- function(silly, IDs){
  
  if(purrr::is_empty(silly)){
    
    stop( paste0("Silly is empty - is this expected?")) # check to see if silly contains anything
    
  }
  
  my.gcl <- get(paste0(silly, ".gcl"), pos = 1)
  
  if(purrr::is_character(IDs)) {IDs <- as.numeric(IDs)}
  
  if(!all(IDs %in% my.gcl$FK_FISH_ID)){
    
    stop(paste0("No IDs were removed. Some of the IDs do not exist for ", silly, "."))
    
  }

  assign(paste(silly, ".gcl", sep = ""),
         my.gcl %>% dplyr::filter(!FK_FISH_ID %in% IDs),
         pos = 1)
  
  message(paste0(length(IDs), " IDs were removed from ", silly, ".gcl"))
  
  if(dim(get(paste0(silly, ".gcl"), pos = 1))[1]==0){
    
    warning(paste0("All IDs were removed from ", silly, ".gcl \n"))
    
  }
  
  return(tibble::tibble(SILLY_CODE = silly, IDs_Removed= list(IDs) %>% purrr::set_names(silly), is_empty = dim(get(paste0(silly, ".gcl"), pos = 1))[1]==0))
  
}