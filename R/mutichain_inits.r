mutichain_inits <- function(npops, nchains, prop = 0.9, type = c("BAYES", "rubias")[1], sillyvec = NULL){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # This function creates either a matrix of initial start values for multiple BAYES chains 
  # or a list of tibbles (length = nchain) containing initial start values for multiple rubias chains.
  # The values for each chain add to 1 and higher start values are given to a different range of pops for each chain.
  # For example, if you have 500 pops and need starting values for 5 chains and prop = 0.9,
  # then the first 100 pops will have starting values = .9/100 and the rest
  # of the pops will have start values = .1/400.  Each chain will have a different 
  # set of pops with this higher start value. 
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #  npops - number of populations
  #
  #  nchains - number of BAYES or rubias chains that you need starting values for.
  #
  #  prop - start value proportion you want to give npops/nchains, default is .9
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
  # 
  # BAYES_inits <- mutichain_inits(npops = length(sillyvec), nchains = 5, prop = 0.9, type = "BAYES", sillyvec = NULL)
  #
  # rubias_inits <- mutichain_inits(npops = length(sillyvec), nchains = 5, prop = 0.9, type = "rubias", sillyvec = sillyvec)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  if(!type%in%c("BAYES", "rubias")){
    
    stop("Only 'BAYES' and 'rubias' types are currently supported")
    
  }
  
  if(nchains < 2){
    
    stop("nchains must be greater than 1")
    
  }
  
  if(type == "rubias" & is.null(sillyvec)){
    
    stop("The user must supply a sillyvec when type = 'rubias'")
    
  }

  r <- npops%%nchains

  avgpop <- floor(npops/nchains)

  groupvec <- c(rep(1:(nchains), each = avgpop), rep(nchains, r)) 

  GroupWeights <- array((1-prop)/(nchains-1), c(nchains, nchains))

  diag(GroupWeights) <- prop  
  
  if(type == "BAYES"){
    
    initmat <- apply(GroupWeights, 1, function(groupweights){
      
      create_prior(groupvec = groupvec, groupweights = groupweights, minval = 0)
      
      })  
    
    dimnames(initmat) <- list(1:npops, paste0("Chain", 1:nchains))
    
  } 
  
  if(type == "rubias"){
    
    initmat <- lapply(1:nchains, function(chain){
      
     tibble::tibble(collection = sillyvec, pi_init = create_prior(groupvec = groupvec, groupweights = GroupWeights[, chain], minval = 0) %>% as.numeric())
        
    }) %>% purrr::set_names(paste0("Chain", 1:nchains))
    
  }
 
  return(initmat)
  
}


