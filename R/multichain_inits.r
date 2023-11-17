#' Generate Initialization Matrix for Multiple Chains
#'
#' This function generates an initialization matrix for multiple chains used in
#' Bayesian data analysis or rubias modeling. The matrix specifies the initial
#' values for each chain to facilitate the convergence of the chains to the
#' posterior distribution.
#'
#' @param npops Integer. The total number of populations to consider.
#' @param nchains Integer. The number of chains to be used. Must be greater than 1.
#' @param prop Numeric (optional). The proportion of the prior distribution to
#'   allocate to each group in the chain initialization matrix. Default is 0.9.
#' @param type Character (optional). The type of initialization to be performed.
#'   Possible values are "BAYES" for [BAYES] or "rubias" for [rubias]
#'   modeling. Default is "BAYES".
#' @param sillyvec Numeric vector (optional). Required only when `type` is set
#'   to rubias". The user must supply a vector containing silly values for each
#'   chain. Default is NULL.
#'
#' @return A matrix or a list of tibbles, depending on the chosen `type`, where
#'   rows represent different populations and columns represent different chains.
#'   For "BAYES" type, the matrix contains the initial prior values for each
#'   population in each chain. For "rubias" type, it returns a list of tibbles,
#'   where each tibble contains a "collection" column with the supplied
#'   `sillyvec` and a "pi_init" column containing the initial prior values for
#'   each population in each chain.
#'
#' @details
#' The function supports two types of initialization:
#'
#' - \code{"BAYES"}: Bayesian initialization where the prior distribution is
#'   divided among the chains proportionally based on the \code{prop} argument.
#'   The matrix returned represents the initial prior values for each population
#'   in each chain.
#'
#' - \code{"rubias"}: Rubias initialization where a \code{sillyvec} is required
#'   to provide initial values for each chain. The prior distribution is divided
#'   among the chains based on the \code{GroupWeights} matrix, which represents
#'   the proportion of prior to be allocated to each group in each chain. The
#'   function returns a list of tibbles, each containing the supplied
#'   \code{sillyvec} in the "collection" column and the initial prior values for
#'   each population in each chain in the "pi_init" column.
#'
#' @examples
#' #Bayesian initialization with 5 populations and 3 chains
#' init_matrix_bayes <- GCLr::multichain_inits(npops = 5, nchains = 3, prop = 0.8)
#'
#' # Rubias initialization with 10 populations, 4 chains, and a supplied sillyvec
#' init_matrix_rubias <- GCLr::multichain_inits(npops = 10, nchains = 4, prop = 0.8, type = "rubias",
#'                                      sillyvec = paste0("Pop", 1:10))
#'
#' @export
multichain_inits <- function(npops, nchains, prop = 0.9, type = c("BAYES", "rubias")[1], sillyvec = NULL){    
  
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
      
      GCLr::create_prior(groupvec = groupvec, groupweights = groupweights, minval = 0)
      
      })  
    
    dimnames(initmat) <- list(1:npops, paste0("Chain", 1:nchains))
    
  } 
  
  if(type == "rubias"){
    
    initmat <- lapply(1:nchains, function(chain){
      
     tibble::tibble(collection = sillyvec, pi_init = GCLr::create_prior(sillyvec = sillyvec, groupvec = groupvec, groupweights = GroupWeights[, chain], minval = 0, type = type) %>% dplyr::pull(pi_param))
        
    }) %>% purrr::set_names(paste0("Chain", 1:nchains))
    
  }
 
  return(initmat)
  
}