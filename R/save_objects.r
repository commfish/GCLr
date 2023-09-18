#' Save R Objects to Disk
#'
#' This function saves R objects with either [base::dput()] or [base::saveRDS()].
#'
#' [base::dput()] saves object code to a text file.
#' [base::saveRDS()] saves objects in rds format, which is a compressed format - saves and reads back into R much faster than a text file.
#'
#' @param objects character vector of objects you wish to save
#' @param path character vector of where you want to save objects
#' @param rds logical; if set to `TRUE` the objects will be saved in rds format,
#'        default is `FALSE` for backwards compatibility.
#' 
#' @return All objects are saved as "object.txt" or "object.rds"
#' 
#' @seealso [base::dput()]
#' @seealso [base::saveRDS()]
#'
#' @examples
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' dir.create(path.expand("~/objects"))
#' 
#' save_objects(objects = paste0(sillyvec, ".gcl"), path = path.expand("~/objects"), rds = TRUE)
#' 
#' save_objects(objects = paste0(sillyvec, ".gcl"), path = path.expand("~/objects"), rds = FALSE)
#' 
#' @export
save_objects <- function(objects, path, rds = FALSE) {
  
  if (!all(objects %in% ls(pos = 1))) {
    stop(paste0(
      "These objects:\n",
      paste(setdiff(objects, ls(pos = 1)), collapse = "\n"),
      "\nare not in your workspace, hoser!!!"
    ))
    
  }
  
  empty <- sapply(objects, function(obj) {
    
    x <- get(obj)
    
    if (!is.null(attr(x = x, which = "problems"))) {
      
      attr(x = x, which = "problems") <- NULL
      
    }  # resolves issue with readr::read_csv not playing nicely with dget if attr-problems exist
    
    if (rds == FALSE) {.
      
      dput(x, paste0(path, "/", obj, ".txt"))
      
    } else {
      
      saveRDS(x, paste0(path, "/", obj, ".rds"))
      
    }
    
  })
  
}