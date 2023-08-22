#' Variance Components Calculation
#'
#' This function calculates the variance components of a three-level hierarchy. This function is intended to help with marker selection (i.e., going from a 96 to 24 SNP baseline).
#' 
#' @param sillyvec A character vector specifying the population silly codes for analysis.
#' @param loci A character vector of locus names to include in the analysis. Can include haploid markers and/or combined loci.
#' @param groupvec A numeric vector specifying the reporting group affiliation for each silly code.
#' @param fstatfile Optional character string specifying the path to an existing FSTAT file (".dat"). If provided, the function verifies compatibility with the input data. (default = NULL)
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)  
#' 
#' @details 
#' Component definitions in the output assume that you are running this function on silly codes that are populations (see Value section). 
#' This function could also be run on silly codes that are collections, but the definitions in the output would be different.
#' 
#' @return A tibble containing the four variance components for each locus:
#'   \itemize{
#'     \item P: For populations (reporting groups in this case).
#'     \item S: For subpopulations within populations (populations within reporting groups).
#'     \item I: For individuals within subpopulations (individuals within populations).
#'     \item G: For alleles within individuals (akin to heterozygosity).
#'    }
#'    
##' 
#' @examples
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' groupvec <- GCLr::ex_baseline %>% 
#'   dplyr::group_by(collection) %>% 
#'   dplyr::summarize(group = unique(repunit)) %>% 
#'   dplyr::mutate(group = factor(group) %>% as.numeric()) %>% 
#'   dplyr::pull(group)
#' 
#' loci <- GCLr::ex_baseline[,-c(1:5)] %>%
#'   names() %>%
#'   gsub(pattern = "*\\.1", x = ., replacement = "") %>%
#'   unique()
#' 
#' GCLr::var_comp(sillyvec = sillyvec, loci = loci, groupvec = groupvec, LocusCtl = GCLr::ex_LocusControl)
#' 
#' @references Weir, B.S. 1995. Genetic Data Analysis. Sinauer Associates, Inc. Sunderland, MA. (pg 184)
#' 
#' @export
var_comp <- function(sillyvec, loci, groupvec, fstatfile = NULL, LocusCtl = LocusControl, ncores = parallel::detectCores()){

  names(sillyvec) <- NULL
  
  nloci <- length(loci)
  
  nalleles <- LocusCtl$nalleles[loci]
  
  ploidy <- LocusCtl$ploidy[loci] 
  
  alleles <- LocusCtl$alleles[loci]
  
  my.gcl <- sapply(sillyvec, function(silly){get(paste(silly, ".gcl", sep = ""), pos = 1)}, simplify = FALSE)
  
  n <- GCLr::silly_n(sillyvec)$n
  
  if(!is.null(fstatfile)) {
    
    dat0 <- hierfstat::read.fstat.data(fstatfile)
    
    if(dim(dat0)[1] != sum(n)) {stop(paste(fstatfile, " does not have the same number of fish as 'sillyvec', verify that the 'fstatfile' corresponds to the fish you want to analyze!"))}

    if(dim(dat0)[2] != (nloci + 1)) {stop(paste(fstatfile, " does not have the same number of loci as 'loci', verify that the 'fstatfile' corresponds to the loci you want to analyze!"))}
    
    if(length(unique(dat0$Pop)) != length(groupvec)) {stop("Number of pops in 'fstatfile' is not equal to length(groupvec)!")}
    
    vars <- names(dat0)
    
    dat <- dat0 %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(Region = rep(groupvec, n)) %>% 
      dplyr::select(Region, dplyr::everything())
      
  } else {
    cat("No 'fstatfile' provided, so creating hierfstat data object\n")

    dat0 <- GCLr::create_hierfstat_data(sillyvec = sillyvec, region = groupvec, pop = seq(length(sillyvec)), loci = loci, ncores = ncores, LocusCtl = LocusCtl) %>% 
      tibble::as_tibble()
    
    dat <- dat0 %>% dplyr::select(-spop) %>% 
      dplyr::rename(Region = region, Pop = pop)
    
  }
  
  G <- max(groupvec)
  
  cat("\nCalculating Variance Components for each locus\n", sep = '')
  
  # Multicore loop
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  `%dopar%` <- foreach::`%dopar%`
  
  VarCompByLocus <- foreach::foreach(locus = loci, .packages = c("tidyverse")) %dopar% {
    
    data <- data.frame(dat[ , 1:2], dat[ , locus])
    
    if(ploidy[locus]==2){
      
      vc <- c(locus, hierfstat::varcomp(data, diploid = TRUE)$overall) %>% 
        setNames(c("locus", "P", "S", "I", "G")) %>% 
        t() %>% 
        tibble::as_tibble()
      
    } else{
      
      vc <- c(locus, hierfstat::varcomp(data, diploid = FALSE)$overall) %>% 
        setNames(c("locus", "P", "S", "I")) %>% 
        t() %>% 
        tibble::as_tibble()
      
    }
     
    if(!sum(is.na(vc))){
      
      vc 
      
    }
    
  } %>% dplyr::bind_rows()
  
  parallel::stopCluster(cl)
  
  return(VarCompByLocus)
  
}