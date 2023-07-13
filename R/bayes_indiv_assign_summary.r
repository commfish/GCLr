#' Summarize BAYES Individual Assignment .CLS Files
#'
#' This function summarizes BAYES individual assignment .CLS files.
#'
#' @param group_names Character vector of group names. Length(group_names) should be equal to max(groupvec).
#' @param groupvec Numeric vector indicating reporting group affiliation for each baseline population. Numbers correspond to group_names.
#' @param mixnames Character vector indicating mixture names (i.e., name of the folder where BAYES output files exist).
#' @param maindir Character vector indicating the directory where the mixture folders are located. **Note: The results for each mixture must be in their own folder with the same name as the mixture.**
#' @param nchains Number of chains to summarize (there should be output files for each chain).
#' @param nreps Numeric vector indicating the number of MCMC reps. Can be a vector of length(mixnames) if there are a different number of reps per mixture.
#' @param burn Proportion of iterations to drop from the beginning of each chain. Default is 0.5, which drops the first 50% of iterations.
#' @param thin Number for the 3rd argument for thin in CreateControlFile.GCL. Recommend thin = 100 for IA.
#'
#' @return A tibble with the following variables:
#'   \describe{
#'     \item{mixname}{The name of the mixture}
#'     \item{id}{The individual id}
#'     \item{group_names}{A variable for each group in group_names containing the proportion of thinned MCMC iterations in which an individual was assigned to a specific reporting group}
#'   }
#'
#' @note
#'   \describe{
#'     \item{Individual ids}{Come from FK_FISH_ID if the mixture .gcl object exists, otherwise they are sequential. Mixture objects can be either new- (tibble) or old-style (list) .gcl objects.}
#'     \item{Dropped individuals}{BAYES will drop mixture individuals from the analysis if they contain an allele that does not exist in the baseline. Check the BAYES summary files for each mixture and remove dropped individuals from the mixture .gcl objects before running this function. Otherwise, the function will produce an error message stating that the number of individuals in the mixture object exceeds the number of individuals in the BAYES output.}
#'   }
#'
#' @examples
#'   \dontrun{
#'   group_names <- c("Crescent", "West", "JCL", "SusYen", "Fish", "KTNE", "Kenai", "Kasilof")
#'   groupvec <- c(1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
#'                4, 3, 6, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8)
#'   mixnames <- c("SUCICD16", "SUCICD16corr", "SUCIND16")
#'   maindir <- "V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Mixture/2016 UCIfisheryMixtures/BAYES/Output"
#'   source("V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Mixture/2016 UCIfisheryMixtures/Output/2016mixturegcls.R") # Get .gcl objects
#'   nchains <- 5
#'   nreps <- 40000
#'   burn <- 0.5
#'   thin <- 100
#'   bayes_indiv_assign_summary(group_names, groupvec, mixnames, maindir, nchains, nreps, burn, thin)
#'   }
#'   
#' @export
bayes_indiv_assign_summary <- function(group_names, groupvec, mixnames, maindir, nchains, nreps, burn = 0.5, thin = 100){
  
  names(mixnames) <- mixnames
  
  if(length(nchains) == 1){
    
    nchains <- rep(nchains, length(mixnames))
    
    names(nchains) <- mixnames
    
  }
  
  if(length(nreps) == 1){
    
    nreps <- rep(nreps, length(mixnames))
    
    names(nreps) <- mixnames
    
  }
  
  myoutputdirs <- sapply(mixnames, function(mixname){
    
    paste(maindir, mixname, sep = "/")
    
  }, simplify = FALSE)
  
  myfilenames <- sapply(mixnames, function(mixname){
    
    paste(myoutputdirs[mixname], "/", mixnames[mixname], "Chain", seq(nchains[mixname]), "CLS.CLS", sep = "")
    
  }, simplify = FALSE)
  
  # Loop through mixtures
  results <- lapply(mixnames, function(mix){
    
    testoutput <- readr::read_table(myfilenames[[mix]][1], col_names = FALSE, col_types = readr::cols(.default = readr::col_double()))
    
    mixsampsize <- dim(testoutput)[1]/floor(nreps[mix]/thin)
    
    chains <- paste0("Chain", seq(nchains[mix]))
    
    skip <- mixsampsize*floor(nreps[mix]*burn/thin)
    
    if(exists(paste(mix, ".gcl", sep = ''))) {
      
      my.gcl <- get(paste(mix, ".gcl", sep = ''))
      
      if(!tibble::is_tibble(my.gcl)) {
        
        indnames <- my.gcl$attributes$FK_FISH_ID
        
      } else{
        
        indnames <- my.gcl$FK_FISH_ID
        
      }
      
    } else {
      
      indnames <- seq(mixsampsize)
      
    }
    
    if(length(indnames) > mixsampsize){
      
      stop(paste0("The number of individuals in ", mix, ".gcl exceeds the number of individuals in the BAYES output. 
                 Check the BAYES summary files to see if individuals were dropped from the analysis because they had 
                 a type or allele not observed in any of the baseline stocks. Remove the dropped individuals from the
                 mixture .gcl objects before running this function."))
    }
    
    # Loop through chains
    output <- lapply(chains, function(chain){
      
      readr::read_table(myfilenames[[mix]][match(chain, chains)], col_names = FALSE, skip = skip, col_types = readr::cols(.default = readr::col_double()) )
      
    }) %>% dplyr::bind_rows()# End chains loop 
    
    nits <- nrow(output)/mixsampsize
    
    C <- ncol(output)
    
    z <- apply(output, 1, which.max)
    
    tapply(X = z, INDEX = rep(seq(mixsampsize), nits), function(zz){
      
      tibble::as_tibble(tabulate(zz, C)/nits)
      
    }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(group = factor(group_names[groupvec], levels = group_names)) %>% 
      dplyr::group_by(group) %>% 
      dplyr::summarize(dplyr::across(tidyselect:::where(is.numeric), sum), .groups = "drop") %>%
      t() %>% 
      janitor::row_to_names(row_number = 1) %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(mixname = mix, id = indnames) %>% 
      dplyr::select(mixname, id, dplyr::all_of(group_names))
    
  }) %>% dplyr::bind_rows() # End mixture loop
  
  return(results)
  
}