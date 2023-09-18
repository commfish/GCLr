#' Save sillys as R objects
#'
#' This function saves sillys as R objects with [base::dput()] or [base::saveRDS()].
#'
#' @param sillyvec character vector of sillys you wish to save
#' @param path character vector of where you want to save sillys
#' @param rds logical; if set to TRUE the sillys will be saved in rds format,
#'        default is FALSE for backwards compatibility.
#'
#' @return All sillyvec are saved as "sillyvec.txt" or "sillyvec.rds", the .gcl is removed.
#'
#' @seealso [GCLr::save_objects()]
#' @seealso [base::dput()]
#' @seealso [base::saveRDS()]
#' 
#' @examples
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' dir.create(path.expand("~/Raw genotypes"))
#' 
#' save_sillys(sillyvec = sillyvec, path = path.expand("~/Raw genotypes"), rds = FALSE)
#' 
#' save_sillys(sillyvec = sillyvec, path = path.expand("~/Raw genotypes"), rds = TRUE)
#'
#' @export
save_sillys <- function(sillyvec, path, rds = FALSE) {

  sillyvec.gcl <- paste0(sillyvec, ".gcl")
  
  if (!all(sillyvec.gcl %in% ls(, pos = 1))) {
    
    stop(paste0("These sillys:\n", paste(setdiff(sillyvec.gcl, ls(, pos = 1)), collapse = "\n"), "\nare not in your workspace, hoser!!!"))
    
  }
  
  empty <- sapply(sillyvec, function(silly) {
    if (rds == FALSE) {
      
      dput(get(paste0(silly, ".gcl")), paste0(path, "/", silly, ".txt"))
      
    } else {
      
      saveRDS(get(paste0(silly, ".gcl")), paste0(path, "/", silly, ".rds"))
      
    }
    
  })
  
}