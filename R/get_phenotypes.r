get_phenotypes <- function(markerset){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function gets all possible combinations of phenotypes for set of loci. 
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   markerset - is a vector of loci (e.g., c("Ots_vatf-251", "Ots_ZR-575"))
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  #   password = "************"
  #
  #   create_locuscontrol(markersuite = "UCI_Chinook_GTSeq_557SNPs", username = "awbarclay", password = password)
  # 
  #   get_phenotypes(markerset = c("Ots_vatf-251", "Ots_ZR-575"))
  #
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function is called on by CombineLoci.gcl and requires a LocusControl object. Run create_locuscontrol prior to this function.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
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