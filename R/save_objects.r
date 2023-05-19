save_objects <- function(objects, path, rds = FALSE) {
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function saves R objects with either `dput` or 'saveRDS'.
  #
  # dput saves object code a text file
  # saveRDS saves objects in rds format, which is a compressed format - saves and reads back into R much faster than a text file.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   objects - character vector of objects you wish to save
  #
  #   path - character vector of where you want to save objects
  #
  #   rds - logical; if set to TRUE the objects will be saved in rds format
  #         default is FALSE for backwards compatablity.
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   All objects are saved as "objects.txt" or "objects.rds"
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # iris <- iris  # bringing into pos = 1
  # save_objects(objects = "iris", path = "Objects", rds = TRUE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (!all(objects %in% ls(pos = 1))) {
    stop(paste0(
      "These objects:\n",
      paste(setdiff(objects, ls(pos = 1)), collapse = "\n"),
      "\nare not in your workspace, hoser!!!"
    ))
    
  }
  
  empty <- sapply(objects, function(obj) {
    if (rds == FALSE) {
      dput(get(obj), paste0(path, "/", obj, ".txt"))
    } else {
      saveRDS(get(obj), paste0(path, "/", obj, ".rds"))
    }
    
  })
  
}