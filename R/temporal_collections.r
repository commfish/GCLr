temporal_collections <- function(sillyvec, region = NULL, min.samps = 50, sep = "."){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #  Get temporal collections from a vector of pooled sillys for hierarchical ANOVA.
  #  Single-silly pops and sillys with sample sizes < "min.samps" will be removed.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of pooled silly codes (AKA pooled pops) without the ".gcl" extention (e.g. sillyvec <- c("KQUART06.KQUART08.KQUART10")). 
  #
  #   region - optional; a tibble containing two numeric variables: 1) pop_no; the position number for each population in sillyvec
  #                                                                 2) region; the region number for each population in sillyvec
  #
  #   min.samples - the minimum number of samples for including a collection 
  #
  #   sep - the separator used when pooling collections into populations; default is "."
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   Returns a tibble with the following variables: silly; non-pooled silly codes, 
  #                                                  pop; numbers indicating the population affiliation for each collection, 
  #                                                  region; numbers indicating the regional affiliation for each collection 
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Coho/Cook Inlet/2019/2019_Cook_Inlet_coho_baseline/2019_Cook_Inlet_coho_baseline_new.Rdata")
  # 
  #   source(paste0(path.expand("~/R/"), "Functions.R"))#GCL functions
  #   
  #   sillyvec <- Final_Pops$collection
  #     
  #   region <- Final_Pops %>% select(pop = order, region = region)
  #     
  #   temporal_collections (sillyvec = sillyvec, region = region, min.samps = 50, sep = ".")  
  # 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  
  if(!all(sillyvec %in% stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))) {  # Do all sillys exist in the environment?
    
    missing_sillys <- setdiff(sillyvec, stringr::str_remove(string = objects(pattern = "\\.gcl", pos = -1, envir = .GlobalEnv), pattern = "\\.gcl"))
    
    stop(paste0("The following sillys are not in your environment:\n", paste0(missing_sillys, collapse = ", ")))
    
  }
  
  tc <- lapply(1:length(sillyvec), function(pop){
    
    silly <- sillyvec[pop]
    
    tcol0 <- base::strsplit(silly, split = paste0("\\", sep))[[1]]
    
    tcol <- silly_n(tcol0) %>% 
      dplyr::filter(n >= min.samps) %>% 
      pull(silly)
    
    if(length(tcol) < 2){
      
      tibble::tibble(pop = pop, silly = NA_character_)
      
    }else(tibble::tibble(pop = pop, silly = tcol))
    
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(!is.na(silly)) %>% 
    dplyr::select(silly, pop)
  
  if(!is.null(region)){
    
    tc <- tc %>% 
      dplyr::left_join(region, by = c("pop"="pop_no")) %>% 
      dplyr::mutate(region = factor(region, levels = unique(region)) %>% as.numeric())
    
  }
  
  return(tc)
  
}

