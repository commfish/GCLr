save_sillys <- function(sillyvec, path, rds = FALSE) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function saves sillys as R objects with `dput`, it is a wrapper for `dput`.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - character vector of sillys you wish to save
  #   path - character vector of where you want to save sillys
  #   rds - logical; if set to TRUE the sillys will be saved in rds format
  #         default is FALSE for backwards compatablity.
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   All sillyvec are saved as "sillyvec.txt" or "sillyvec.rds", the .gcl is removed
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # iris <- iris  # bringing into pos = 1
  # save_sillys(sillyvec = "SGILL15D15", path = "Raw genotypes")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  sillyvec.gcl <- paste0(sillyvec, ".gcl")
  
  if (!all(sillyvec.gcl %in% ls(, pos = 1))) {
    stop(paste0(
      "These sillys:\n",
      paste(setdiff(sillyvec.gcl, ls(, pos = 1)), collapse = "\n"),
      "\nare not in your workspace, hoser!!!"
    ))
  }
  
  empty <- sapply(sillyvec, function(silly) {
    if (rds == FALSE) {
      dput(get(paste0(silly, ".gcl")), paste0(path, "/", silly, ".txt"))
    } else {
      saveRDS(get(paste0(silly, ".gcl")), paste0(path, "/", silly, ".rds"))
    }
    
  })
  
}