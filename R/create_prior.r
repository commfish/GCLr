#' Create Prior Objects for BAYES or rubias
#'
#' This function creates a prior object for BAYES or rubias.
#'
#' @param groupvec Numeric vector indicating the group affiliation of each population in `sillyvec`.
#' @param groupweights The weights for each group in `groupvec`. For example, to give the same weight to all groups, use \code{groupweights = rep(1/max(groupvec), max(groupvec))}.
#' @param minval the minimum value of the weight given to each population (default = 0.01).
#' @param type Whether you want a prior object for "BAYES" or \pkg{rubias}. Default is "BAYES".
#' @param sillyvec A vector of silly codes without the ".gcl" extension. This argument is needed when \code{type = "rubias"}.
#'
#' @returns  For \code{type = "BAYES"}, the function outputs a vector with names = \code{groupvec}.
#'   For \code{type = "rubias"}, the function outputs a tibble with "collection" and "pi_param" variables.
#'
#' @note The BAYES prior object is supplied to the prior argument in \code{create_bayes_ctl()}.
#'   The rubias prior object is supplied to the pi_prior argument in \code{run_rubias_mix()}.
#'
#' @examples
#' \dontrun{
#' groupvec <- c(1, 1, 2, 2, 3, 4, 4, 5, 2, 5)
#' sillyvec <- c("KQUART060809", "KCRESC06", "KKENU09", "KJUNE050607", "KRUSSR05060708", "KBENJ0506", "KKILL0506", "KFUNN0506", "KKENAI030406", "KSLIK040508")
#' BAYES_prior <- create_prior(groupvec = groupvec, groupweights = rep(1/max(groupvec), max(groupvec)), minval = 0.01, type = "BAYES", sillyvec = NULL)
#' rubias_prior <- create_prior(groupvec = groupvec, groupweights = rep(1/max(groupvec), max(groupvec)), minval = 0.01, type = "rubias", sillyvec = sillyvec)
#' }
#'
#' @details This function is also used by [GCLr::multichain_inits()] to create initial starting values.
#'
#' @export
create_prior <- function(groupvec, groupweights, minval = 0.01, type = c("BAYES", "rubias")[1], sillyvec = NULL){

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