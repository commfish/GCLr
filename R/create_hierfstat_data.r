create_hierfstat_data <- function(sillyvec, region = NULL, pop,loci, ncores  = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   Create a hierfstat data object.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   region - optional; a numeric vector of numbers indicating the regional affiliation for each silly in sillyvec; 
  #                      include this argument when computing hierarchical F statistics.
  #
  #   pop - a numeric vector indicating the pop affiliation for each silly in sillyvec.
  #         if there is only one silly per population in in sillyvec, then this is a sequence of numbers the length of sillyvec - i.e., seq(length(sillyvec)).
  #
  #   loci - a character vector of locus names
  #
  #   ncores - a numeric vector of length one indicating the number of cores to use
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   Returns a hierfstat data object containing region (if supplied), population, and sub population numbers and genotypes in single-column format.
  #   
  #   Note: leading zeros are dropped from the genotypes (e.g., 0102 == 102), but hierfstat functions account for this.
  #
  #   Note: If is.null(region), the object will not contain a region column.  The region and sub population columns are intended for computing hierarchical F statistics (i.e., hierfstat::varcomp, hierfstat::varcomp.glob), 
  #         and they will need to be removed from the object if supplying to hierfstat functions that compute population-level statistics.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Coho/Cook Inlet/2019/2019_Cook_Inlet_coho_baseline/2019_Cook_Inlet_coho_baseline_new.Rdata")
  # 
  #   source(paste0(path.expand("~/R/"), "Functions.R"))#GCL functions
  # 
  #   sillyvec104 <- Final_Pops$collection
  # 
  #   region <- Final_Pops %>% select(pop = order, region = region)
  # 
  #   temporal_collections <- temporal_collections(sillyvec = sillyvec104, region = region, min.samps = 50, sep = ".")
  # 
  #   fstat.dat <- create_hierfstat_data(sillyvec = temporal_collections$silly, region = temporal_collections$region, pop = temporal_collections$pop, loci = loci81, ncores = 8)
  # 
  #   hierfstat::varcomp.glob(levels = fstat.dat[,1:3], loci = fstat.dat[,-c(1:3)])
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(!exists("LocusControl")){
    
    stop("'LocusControl' not yet built.")
    
  }
  
  if(sum(is.na(match(loci,LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci,LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  if(!all(sillyvec %in% stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))) {  # Do all sillys exist in the environment?
    
    missing_sillys <- setdiff(sillyvec, stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))
    
    stop(paste0("The following sillys are not in your environment:\n", paste0(missing_sillys, collapse = ", ")))
    
  }
  
  alleles <- LocusControl$alleles[loci] %>% 
    dplyr::bind_rows(.id = "locus")
  
  my.gcl <- lapply(sillyvec, function(silly){
    
    get(paste(silly,".gcl", sep = ""), pos = 1)
    
    }) %>% 
    purrr::set_names(sillyvec)
  
  # Multicore loop
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  output <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse")) %dopar% {
    
    new.gcl <- my.gcl[[silly]]
    
    spop_no <- match(silly, sillyvec)
    
    pop_no <- pop[spop_no]
    
    if(is.null(region)){
      
      region_no <- NA_integer_
      
    }else{region_no <- region[spop_no]}
    
    scores_names <- sapply(loci, function(locus) {
      
      c(locus, paste0(locus, ".1"))
      
      }) %>% as.vector() 
    
    scores <- new.gcl %>% 
      dplyr::select(tidyselect::all_of(scores_names)) %>% 
      as.data.frame(stringsAsFactors = FALSE)
    
    lapply(loci, function(loc){
      
      variables <- c(loc, paste(loc, 1, sep = "."))
      
      my.alleles <- alleles %>% 
        dplyr::filter(locus==loc)
      
      scores %>%
        dplyr::select(tidyselect::all_of(variables)) %>% 
        dplyr::mutate(across(dplyr::everything(), .fns = ~factor(., levels = my.alleles$call))) %>% 
        dplyr::mutate(across(dplyr::everything(), .fns = ~as.numeric(.))) %>%
        dplyr::mutate(across(dplyr::everything(), .fns = ~stringr::str_pad(., width = 2, side = "left", pad = "0"))) %>% 
        tidyr::unite(col = !!rlang::as_name(loc), tidyselect::all_of(variables), sep = "", na.rm = TRUE) %>% 
        dplyr::mutate(across(dplyr::everything(), .fns = ~as.numeric(.)))
      
      
    }) %>% 
      dplyr::bind_cols() %>% 
      purrr::set_names(loci) %>% 
      dplyr::mutate(region = region_no, pop = pop_no, spop = spop_no) 
  } %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(`region`, pop, spop, dplyr::everything())# end multicore loop
  
  parallel::stopCluster(cl)
  
  if(is.null(region)) {
    
    output <- output %>%
      dplyr::select(-`region`)
    
  }
  
  return(output)
  
}  
