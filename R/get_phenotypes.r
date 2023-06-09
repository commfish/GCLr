#' Get All Combinations of Phenotypes
#' 
#' This function gets all possible combinations of phenotypes for set of loci.
#' 
#' @param markerset is a vector of loci (e.g., c("Ots_vatf-251", "Ots_ZR-575"))
#' 
#' @param LocusControl an object created by [GCLr::create_locuscontrol()]
#' 
#' @note this function is called on by [GCLr::combine_loci()] and requires a LocusControl object. Run [GCLr::create_locuscontrol()] prior to this function.
#' 
#' @return a vector of phenotypes
#'
#' @examples
#' 
#' GCLr::get_phenotypes(markerset = c("One_vatf-214", "One_ZNF-61"), LocusControl = GCLr::ex_LocusControl)
#' 
#' @export
get_phenotypes <- function(markerset, LocusControl = LocusControl){
  
  alleles <- LocusControl$alleles[markerset]
  
  ploidy <- LocusControl$ploidy[markerset]
  
  newgntps <- lapply(
    
    lapply(markerset, function(marker){
      
      t(combn(rep(alleles[[marker]]$call, ploidy[marker]), ploidy[marker]))
      
    }), function(gmat){
      
      na <- NULL
      
      nr <- nrow(gmat)
      
      nc <- ncol(gmat)
      
      for(i in 1:nr){
        
        if(nc > 1){
          
          vec <- sort(gmat[i, ])
          
        } else{
          
          vec <- gmat[i]
          
        } 
        
        na <- c(na, paste(vec, collapse = ""))
        
      }
      
      sort(unique(na))
      
    })
  
  if(length(unique(markerset)) == 1){
    
    newalleles = newgntps[[1]]
    
  } else {
    
    temp <- expand.grid(newgntps, stringsAsFactors = FALSE)     
    
    nr <- nrow(temp)
    
    na <- NULL
    
    for(i in 1:nr){
      
      na <- c(na, paste(as.character(temp[i, ]), collapse = ""))
      
    }      
    
    newalleles <- sort(unique(na))
    
  }  # else 
  
  return(newalleles)
}