#' @title Create Pairwise F~CT~ Matrix
#'
#' @description 
#' This function generates a matrix of pairwise Fct values.
#'
#' @param sillyvec A character vector of silly codes without the ".gcl" extension.
#' @param loci A character vector of locus names.
#' @param group_names character vector of group names the length of \code{max(groupvec)}
#' @param groupvec a numeric vector indicating the group affiliation of each pop in `sillyvec`
#' @param ncores A numeric value for the number of cores to use in a \pkg{foreach} `%dopar%` loop (default = 4). 
#' The number of cores cannot exceeds the number on your device ([parallel::detectCores()]).
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()] (default = LocusControl).
#'
#' @returns A matrix of F~CT~ values between reporting groups
#'
#' @details
#' The function calls on [GCLr::create_hierfstat_data()] to create a \pkg{hierfstat} data object internally and the uses[hierfstat::varcomp()] to perform a 4-level analysis of variance (NOVA).
#' The variance components are then used to calculate F~CT~ between each reporting group. F~CT~ is a measure of genetic differentiation between groups of populations and is analogous to Fst between populations. 
#'
#' @examples
#'   
#'   sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#'   group_names <- GCLr::ex_baseline %>%
#'     dplyr::group_by(collection) %>%
#'     dplyr::summarize(group = unique(repunit)) %>%
#'     dplyr::pull(group) %>%
#'     unique()
#'   
#'   groupvec <- GCLr::ex_baseline %>%
#'     dplyr::group_by(collection) %>%
#'     dplyr::summarize(group = unique(repunit)) %>%
#'     dplyr::mutate(group = factor(group, levels = group_names) %>% as.numeric()) %>%
#'     dplyr::mutate(collection = factor(collection, levels = sillyvec)) %>% 
#'     dplyr::arrange(collection) %>% 
#'     dplyr::pull(group)
#' 
#' 
#'   loci <- GCLr::ex_baseline[,-c(1:5)] %>%
#'     names() %>%
#'     gsub(pattern = "*\\.1", x = ., replacement = "") %>%
#'     unique()
#'   
#'   
#'  pairwise_fct(sillyvec = sillyvec, loci = loci, group_names = group_names,  groupvec = groupvec, ncores = 20, LocusCtl = GCLr::ex_LocusControl)
#' 
#' @export
pairwise_fct <- function(sillyvec, loci, group_names, groupvec, ncores = parallel::detectCores(), LocusCtl = LocusControl) {
    
    if(ncores > parallel::detectCores()) {
      
      stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
      
    }
    
    if (.Platform$OS.type == "windows"){utils::flush.console()}
    
    start.time <- Sys.time() 
    
    n <- GCLr::silly_n(sillyvec) %>% 
      dplyr::pull(n) %>% 
      purrr::set_names(sillyvec)
    
    nsillys <- length(sillyvec)
    
    nloci <- length(loci)
    
    ngroups <- length(group_names)
    
    message("\nCreate hierfstat data object\n", sep = '')
    
    dat <- GCLr::create_hierfstat_data(sillyvec = sillyvec, region = groupvec, pop = seq_along(sillyvec), loci = loci, ncores = ncores, LocusCtl = LocusCtl)
    
    pairs <- combn(group_names, 2)
    
    pairnames <- apply(pairs, 2, function(grp){
      
      paste(grp, collapse = ".")
      
    })
    
    dimnames(pairs)[[2]] <- pairnames
    
    ploidy <- LocusCtl$ploidy
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Calculate variance components for each pair of sillys
    message("\nCalculate variance components for each pair of groups\n", sep = '')
    
    cl <- parallel::makePSOCKcluster(ncores)
    
    doParallel::registerDoParallel(cl, cores = ncores)  
    
    `%dopar%` <- foreach::`%dopar%`
    
    vc <- foreach::foreach(pair = pairnames, .packages = c("hierfstat", "magrittr")) %dopar% {
      
      grps <- pairs[, pair]
      
      grpIND <- sapply(group_names, function(grp){match(grp, group_names)})
      
      mydat <- dat %>% 
        dplyr::filter(region %in% grpIND[grps]) %>% 
        dplyr::select(region, pop, dplyr::all_of(loci))
      
      mito.loci <- which(ploidy[loci] == 1)
      
      # Create matrix of variance components with order loci
      if(length(mito.loci) == 0) {
        
        # Diploid loci
        vc.pair <- t(sapply(loci, function(locus) {
          
          mydata <- mydat %>% 
            dplyr::select(region, pop, dplyr::all_of(locus))
          
          mydata_no_na <- dplyr::filter(.data = mydata, !is.na(!!dplyr::sym(locus)))
          
          if (dplyr::n_distinct(mydata_no_na$region) == 2) {
            
            hierfstat::varcomp(data = mydata, diploid = TRUE)$overall
            
          } else {
            
            c(NA, NA, NA, NA)
            
          }
          
        }))
        
      } else {
        
        # Diploid loci
        vc.pair <- rbind(t(sapply(loci, function(locus) {
          
          mydata <- mydat %>% 
            dplyr::select(region, pop, dplyr::all_of(locus))
          
          mydata_no_na <- dplyr::filter(.data = mydata, !is.na(!!dplyr::sym(locus)))
          
          if (dplyr::n_distinct(mydata_no_na$region) == 2) {
            
            hierfstat::varcomp(data = mydata, diploid = TRUE)$overall
            
          } else {
            
            c(NA, NA, NA, NA)
            
          }
          
        })),
        
        # Haploid loci
        t(sapply(loci[mito.loci], function(locus) {
          
          mydata = mydat %>% dplyr::select(region, pop, dplyr::all_of(locus))
          
          mydata_no_na <- dplyr::filter(.data = mydata, !is.na(!!dplyr::sym(locus)))
          
          if (dplyr::n_distinct(mydata_no_na$levels) == 2) {
            
            c(hierfstat::varcomp(data = data.matrix(mydata), diploid = FALSE)$overall,0)
            
          } else {
            
            c(NA, NA, NA, NA)
            
          }
          
        } ))
        
        )[loci, ]
        
      }
      
      # Replace NAs with 0s
      vc.pair[is.na(vc.pair)] <- 0
      vc.pair
      
    } # End VC calc
    
    names(vc) <- pairnames
    
    # Calculate Fct from VC for each pair
    Fct.ls <- lapply(vc, function(pair) {
      
      sum(pair[, 1]) / sum(pair)
      
    })
    
    # Create matrix of pairwise Fst
    Fct <- array(0, c(length(group_names), length(group_names)), dimnames = list(group_names, group_names))
    
    Fct[lower.tri(Fct)] <- unlist(Fct.ls)
    
    Fct[upper.tri(Fct)] <- t(Fct)[upper.tri(Fct)]
    
    # Stop cluster
    parallel::stopCluster(cl)
    
    stop.time <- Sys.time()
    
    fulltime <- stop.time - start.time
    
    print(fulltime) 
    
    return(Fct)
    
  }