random_inits <- function(groupvec, groupweights = rep(1/max(groupvec), max(groupvec)), nchains = 5,  type = c("BAYES", "rubias")[1], sillyvec = NULL){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # This function creates either a matrix of randomly generated initial start values for multiple BAYES chains 
  # or a list of tibbles (length = nchain) containing randomly generated initial start values for multiple rubias chains.
  # The values for each chain add to 1.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #  groupvec - numeric vector indicating the group affiliation of each pop in sillyvec
  #
  #  groupweights - the weights for each group in groupvec. For example, to give the same weight to all groups groupweights = rep(1/max(groupvec), max(groupvec))
  #
  #  nchains - number of BAYES or rubias chains that you need starting values for.
  #
  #  type - whether you want an initial start value object for "BAYES" or "rubias"
  #
  #  sillyvec - a vector of silly codes without the ".gcl" extension, this argument is needed when type = "rubias"
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  # For type = "BAYES" the function outputs a matrix with nrow = npops and ncol = nchains.
  # For type = "rubias" the function outputs a list of length = nchains, and each element is a tibble
  # with "collection" and "pi_init" variables. 
  #
  # Note: The BAYES initial start value object is supplied to the initmat argument in create_bayes_ctl()
  #       The rubias initial start value object is supplied to the pi_init argument in run_rubias_mix(); however,[ 
  #       the user must subset the list to supply a single tibble for each rubias chain. 
  #       For example, use rubias_inits$Chain1, rubias_init[["Chain1"]] or rubias_init[[1]] to subset for chain 1.
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # sillyvec <- c("KQUART060809", "KCRESC06", "KKENU09", "KJUNE050607", "KRUSSR05060708", "KBENJ0506", "KKILL0506", "KFUNN0506", "KKENAI030406", "KSLIK040508")
  # groupvec <- c(1, 1, 2, 2, 3, 4, 4, 5, 2, 5)
  # 
  # BAYES_inits <- random_inits(groupvec = groupvec, groupweights = rep(1/max(groupvec), max(groupvec)), nchains = 5, type = "BAYES", sillyvec = NULL)
  #
  # rubias_inits <- random_inits(groupvec = groupvec, groupweights = rep(1/max(groupvec), max(groupvec)), nchains = 5, type = "rubias", sillyvec = sillyvec)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  if(!type %in% c("BAYES", "rubias")){
    
    stop("Only 'BAYES' and 'rubias' types are currently supported")
    
  }
  
  if(nchains < 2){
    
    stop("nchains must be greater than 1")
    
  }
  
  if(type == "rubias" & is.null(sillyvec)){
    
    stop("The user must supply a sillyvec when type = 'rubias'")
    
  }
  
  chains <- paste0("Chain", 1:nchains)

  Cg <- table(groupvec)
  
  popweights <- groupweights[groupvec]/Cg[groupvec]

  if(type == "BAYES"){
    
    ans <- sapply(chains, function(chain){
      
      sapply(popweights, function(popweight){
        
        rgamma(1,popweight, 1)
        
        })
      
      })
    
    ans <- apply(ans, 2, function(x){x/sum(x)})
    
  }
  
  if(type == "rubias"){
    
    ans <- lapply(chains, function(chain){
      
     x <- sapply(popweights, function(popweight){
        
        rgamma(1, popweight, 1)
        
      })
     
     tibble::tibble(collection = sillyvec, pi_init = x/sum(x))
     
    }) %>% purrr::set_names(paste0("Chain", seq(nchains)))
 
  }

  return(ans)
  
}

