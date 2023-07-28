#' Save R objects to disk
#'
#' This function saves R objects with either `dput` or `saveRDS`.
#'
#' `dput` saves object code to a text file.
#' `saveRDS` saves objects in rds format, which is a compressed format - saves and reads back into R much faster than a text file.
#'
#' @param objects character vector of objects you wish to save
#' @param path character vector of where you want to save objects
#' @param rds logical; if set to TRUE the objects will be saved in rds format,
#'        default is FALSE for backwards compatibility.
#' 
#' @return All objects are saved as "objects.txt" or "objects.rds"
#'
#' @examples
#' iris <- iris  # bringing into pos = 1
#' save_objects(objects = "iris", path = "Objects", rds = TRUE)
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
    if (rds == FALSE) {
      dput(x, paste0(path, "/", obj, ".txt"))
    } else {
      saveRDS(x, paste0(path, "/", obj, ".rds"))
    }
    
  })
  
}