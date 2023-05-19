create_prior <- function(groupvec, groupweights, minval = 0.01, type = c("BAYES", "rubias")[1], sillyvec = NULL){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # This function creates a prior objects for BAYES or rubias.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #  groupvec - numeric vector indicating the group affiliation of each pop in sillyvec
  #
  #  groupweights - the weights for each group in groupvec. For example, to give the same weight to all groups groupweights = rep(1/max(groupvec), max(groupvec))
  #
  #  minval - the minimum value of the weight given to each population. 
  #
  #  type - whether you want a prior object for "BAYES" or "rubias"
  #
  #  sillyvec - a vector of silly codes without the ".gcl" extension, this argument is needed when type = "rubias"
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  # For type = "BAYES" the function outputs a vector with names = groupvec.
  # For type = "rubias" the function outputs a tibble with "collection" and "pi_param" variables. 
  #
  # Note: The BAYES prior object is supplied to the prior argument in create_bayes_ctl()
  #       The rubias prior object is supplied to the pi_prior argument in run_rubias_mix()
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # groupvec <- c(1, 1, 2, 2, 3, 4, 4, 5, 2, 5)
  #
  # sillyvec <- c("KQUART060809", "KCRESC06", "KKENU09", "KJUNE050607", "KRUSSR05060708", "KBENJ0506", "KKILL0506", "KFUNN0506", "KKENAI030406", "KSLIK040508")
  #
  # BAYES_prior <- create_prior(groupvec = groupvec, groupweights = rep(1/max(groupvec), max(groupvec)), minval = 0.01, type = "BAYES", sillyvec = NULL)
  #
  # rubias_prior <- create_prior(groupvec = groupvec, groupweights = rep(1/max(groupvec), max(groupvec)), minval = 0.01, type = "rubias", sillyvec = sillyvec)
  # 
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # This function is also used by MulitChainInits.GCL() to create initial starting values.
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!type%in%c("BAYES", "rubias")){
    
   stop("Only 'BAYES' and 'rubias' types are currently supported")
    
  }
  
  if(type == "rubias" & is.null(sillyvec)){
    
    stop("The user must supply a sillyvec when type = 'rubias'")
    
  }
  
  Cg <- table(groupvec)

  IND <- groupweights < minval

  groupweights[!IND] <- (1-minval*sum(IND))*groupweights[!IND]/sum(groupweights[!IND])  

  groupweights[IND] <- minval

  if(type == "BAYES"){
    
    popweights <- groupweights[groupvec]/Cg[groupvec]
  
  }
  
  if(type == "rubias"){
    
    popweights <- tibble::tibble(collection = sillyvec, pi_param = as.numeric(groupweights[groupvec]/Cg[groupvec]))
    
  }

  return(popweights)
  
}