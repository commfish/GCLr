dupcheck_within_silly <- function(sillyvec, loci = LocusControl$locusnames, quantile = NULL, minnonmissing = 0.6, minproportion = 0.95, ncores = 4, LocusCtl = LocusControl){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function checks for duplicate individuals within each silly in "sillyvec".
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   loci - vector of locus names; if set to NULL all loci in the ".gcl" objects will be used.
  #
  #   quantile - this argument along with minproportion are used together to determine the cut-off proportion at which a pair of duplicates 
  #                                is defined: i.e. proportion = max(quantile(duplication, quantile), minproportion. 
  #                                Setting "quantile" equal to NULL will skip the calculation of the duplication distribution and will run much faster
  #
  #   minnonmissing - the proportion of loci that a pair must share non missing in order to be reported
  #
  #   minproportion - the proportion of shared non-missing loci that must be shared between the indivdiuals to be reported as a matching pair, this is passed through to rubias as well 
  #
  #   ncores - the number of cores to use in a foreach %dopar% loop. If the number of core exceeds the number on your device, then ncores defaults to detectCores()
  # 
  #   LocusCtl - an object created by [GCLr::create_locuscontrol()], (default = LocusControl) 
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    When quantile is set to NULL, returns a tibble of duplicate pairs of individuals by silly.
  #    When quantile is a number, a list containing a tibble of duplicate pairs of individuals by silly and tibble the porportion of duplication for each pair of individuals.
  #
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # password <- "************"
  # username <- "******"
  # create_locuscontrol(markersuite = "Sockeye2011_96SNPs", username = username, password = password)
  # sillyvec <- c("SMCDO03", "SNEVA13")
  # loki2r(sillyvec = sillyvec, username = username, password = password)
  # remove_ind_miss_loci(sillyvec = sillyvec)
  # pool_collections(c("SMCDO03", "SNEVA13"))
  # 
  # dupcheck <- dupcheck_within_silly(sillyvec = "SMCDO03.SNEVA13", loci = LocusControl$locusnames, quantile = 0.99, minproportion = 0.95, ncores = 4)
  # dupcheckNULLQantile <- dupcheck_within_silly(sillyvec = "SMCDO03.SNEVA13", loci = LocusControl$locusnames, quantile = NULL, minnonmissing = 0.6, minproportion = 0.95, ncores = 4)
  # 
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   When quantile is set to NULL this function utilizes rubias::close_matching_samples() to perform the duplicate check and it much faster than when you set a quantile.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  start.time <- Sys.time() 
  
  if(!all(loci %in% LocusCtl$locusnames)){
    
    stop(paste0("The following `loci` were not found in `LocusControl`:\n", paste(setdiff(loci, LocusCtl$locusnames), collapse = "\n")))
    
  }
  
  if(ncores > parallel::detectCores()) {
    
    ncores <- parallel::detectCores()  
    
    }
  
  nsilly <- length(sillyvec)
  
  nloci <- length(loci)
  
  ploidy <- LocusCtl$ploidy[loci]
  
  scores_cols <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
    as.vector()  # This keeps the scores columns in the correct order when there are loci with similar names.
   
  attr <-
    c(
      "FK_FISH_ID",
      "COLLECTION_ID",
      "SILLY_CODE",
      "PLATE_ID",
      "PK_TISSUE_TYPE",
      "CAPTURE_LOCATION",
      "CAPTURE_DATE",
      "END_CAPTURE_DATE",
      "MESH_SIZE",
      "MESH_SIZE_COMMENT",
      "LATITUDE",
      "LONGITUDE",
      "AGENCY",
      "VIAL_BARCODE",
      "DNA_TRAY_CODE",
      "DNA_TRAY_WELL_CODE",
      "DNA_TRAY_WELL_POS",
      "CONTAINER_ARRAY_TYPE_ID",
      "SillySource"
    )
  
  sillyvec_new <- silly_n(sillyvec) %>% 
    dplyr::filter(n > 1) %>% 
    dplyr::pull(silly)  # Removing single individual collections from sillyvec.
  
  my.gcl <- sapply(sillyvec_new, function(silly){
    
    gcl <- get(paste0(silly, ".gcl"), pos = 1) %>%
      dplyr::select(tidyselect::all_of(attr),
                    tidyselect::all_of(scores_cols))
    
    maxna <- max(rowSums(is.na(gcl[, scores_cols])))  # What is max number of missing loci?
    
    if(maxna == 2 * nloci){
      
      stop(paste0("Some of the individuals in sillyvec have no data for all loci. Run remove_ind_miss_loci(sillyvec, loci) prior to checking for duplicate individuals."))
      
    }
    
    # Added this if statement code for haploid markers, rubias::close_matching_samples() was counting them as missing loci because of the NAs in the allele2 column. 
    # This can be removed Eric Anderson fixes the function. 
    if(any(ploidy == 1)) {
      
      haploci <- names(ploidy[ploidy == 1])
      
      for(locus in haploci){
        
        gcl[[paste0(locus, ".1")]] <- gcl[[locus]]
        
      }
      
      gcl 
      
    }
    
    gcl
    
  }, simplify = FALSE) 
  
  `%dopar%` <- foreach::`%dopar%`
  
  # Start if NULL quantile
  if(is.null(quantile)){
    
    ## Loop through sillys
    cl <- parallel::makePSOCKcluster(ncores)
    
    doParallel::registerDoParallel(cl, cores = ncores)  
    
    dupcheck0 <- foreach::foreach(silly = sillyvec_new, .export = c("loci"), .packages = c("tidyverse","rubias")) %dopar% {
      
      new.gcl <- my.gcl[[silly]] %>%
        dplyr::mutate(
          sample_type = "reference",
          repunit = NA_character_,
          collection = SILLY_CODE
        ) %>%
        tidyr::unite(col = "indiv", SILLY_CODE, FK_FISH_ID, sep = "_") %>% 
        dplyr::select(sample_type,
                      repunit,
                      collection,
                      indiv,
                      tidyselect::all_of(scores_cols))
      
      rubias::close_matching_samples(D = new.gcl, gen_start_col = 5, min_frac_non_miss = minnonmissing, min_frac_matching = minproportion)
      
    } %>% dplyr::bind_rows()  # End multicore loop
    
    parallel::stopCluster(cl)
    
    dupcheck <- dupcheck0 %>% 
      tidyr::separate(indiv_1, into = c(NA, "ID1"), sep = "\\_(?=[^\\_]+$)", extra = "drop") %>% 
      tidyr::separate(indiv_2, into = c(NA, "ID2"), sep = "\\_(?=[^\\_]+$)", extra = "drop") %>%
      dplyr::mutate(silly = collection_1, 
                    proportion = num_match/num_non_miss) %>% 
      dplyr::select(silly, ID1, ID2, proportion)
    
    # Calculate the proportion of missing scores for ID1 and ID2 and create the output report.
    
    report <- lapply(unique(dupcheck$silly), function(silly){
      
      new.gcl <- my.gcl[[silly]]
      
      dups <- dupcheck %>% 
        dplyr::filter(silly==!!silly) %>% 
        dplyr::mutate(order = seq(length(silly)), 
                      ID1 = as.numeric(ID1), 
                      ID2 = as.numeric(ID2)
        )
      
      ID1 <- new.gcl %>% 
        dplyr::mutate(ID1 = as.numeric(FK_FISH_ID)) %>% 
        dplyr::filter(ID1 %in% dups$ID1) %>% 
        dplyr::select(ID1, tidyselect::all_of(loci)) %>% 
        tidyr::gather(-ID1, key = "Locus", value = "allele") %>% 
        dplyr::group_by(ID1) %>% 
        dplyr::summarize(Missing1 = sum(is.na(allele)), .groups = "drop_last") %>% 
        dplyr::left_join(dups, by = "ID1")%>% 
        dplyr::arrange(order)
      
      ID2 <- new.gcl %>% 
        dplyr::mutate(ID2 = as.numeric(FK_FISH_ID)) %>% 
        dplyr::filter(ID2 %in% dups$ID2) %>% 
        dplyr::select(ID2, tidyselect::all_of(loci)) %>% 
        tidyr::gather(-ID2, key = "Locus", value = "allele") %>% 
        dplyr::group_by(ID2) %>% 
        dplyr::summarize(Missing2 = sum(is.na(allele)), .groups = "drop_last") %>% 
        dplyr::left_join(dups, by = "ID2") %>% 
        dplyr::arrange(order)
      
      dplyr::full_join(ID1, ID2, by = c("ID1", "silly", "ID2", "proportion", "order")) %>% 
        dplyr::select(silly, ID1, ID2, Missing1, Missing2, proportion)
      
    }) %>%  dplyr::bind_rows()  
    
    print(Sys.time() - start.time)
    
    return(report)
    
  }  # End if NULL quantile
  
  # Start if quantile 
  if(!is.null(quantile)){
    
    cl <- parallel::makePSOCKcluster(ncores)
    
    doParallel::registerDoParallel(cl, cores = ncores)  
    
    # multicore loop
    dupcheck <- foreach::foreach(silly = sillyvec, .packages = "tidyverse") %dopar% {
      
      new.gcl <- my.gcl[[silly]]
      
      IDs <- new.gcl$FK_FISH_ID
      
      n <- dim(new.gcl)[1]
      
      nloci <- length(loci)
      
      # Combine allele 1 and 2 into single columns for each locus separated by a period.
      scores0 <- new.gcl[ , c("FK_FISH_ID", scores_cols)] 
      
      if(n < 2){
        
        report <- tibble::tibble(silly = silly, 
                                 ID1 = NA, 
                                 ID2 = NA, 
                                 Missing1 = NA, 
                                 Missing2 = NA, 
                                 proportion = NA)  # If only one individual the results will be a tibble of NAs. Can't do next() using applys.
        
      } else {
        
        ncombs <- choose(n, 2)
        
        mycombs <- combn(IDs, 2)
        
        scores0_allele_1 <- scores0 %>% 
          dplyr::select(FK_FISH_ID, !!!loci) %>% 
          tidyr::gather(locus, allele_1, -FK_FISH_ID)
        
        scores0_allele_2 <- scores0 %>% 
          dplyr::select(FK_FISH_ID, !!!paste0(loci, ".1")) %>% 
          tidyr::gather(locus, allele_2, -FK_FISH_ID) %>% 
          dplyr::mutate(locus = stringr::str_remove(string = locus, pattern = "\\.1"))
        
        scores1 <- dplyr::left_join(x = scores0_allele_1, y = scores0_allele_2, by = c("FK_FISH_ID", "locus")) %>% 
          tidyr::unite(col = "genotype", c("allele_1", "allele_2"), sep = ".") %>% 
          dplyr::mutate(genotype = dplyr::na_if(genotype, "NA.NA")) %>% 
          tidyr::spread(locus, genotype) %>% 
          tibble::column_to_rownames("FK_FISH_ID")
        
        duplication <- sapply(1:ncombs, function(comb){
          
          compair <- scores1[as.character(mycombs[1, comb]), ] == scores1[as.character(mycombs[2, comb]), ]
          
          sum(compair[!is.na(compair)]) / sum(!is.na(compair))
          
        }) %>% purrr::set_names(sapply(1:ncombs, function(comb){
          
          paste(mycombs[,comb], collapse=".")
          
        }))
        
        proportion <- max(quantile(duplication, quantile), minproportion)
        
        dupIND <- duplication >= proportion  # Made this greater than or equal to so it matches the output from rubias::close_matching_samples()
        
        if(sum(dupIND)){
          
          dups <- data.frame(ID1 = mycombs[1, dupIND], ID2 = mycombs[2, dupIND])
          
          missing <- t(sapply(1:nrow(dups), function(dup){
            
            vec <- match(as_vector(dups[dup, ]), IDs)
            
            c("Missing1" = sum(is.na(scores1[vec[1], ])),
              "Missing2" = sum(is.na(scores1[vec[2], ])))
            
          }, simplify = TRUE))
          
          report <- dplyr::bind_cols(dplyr::as_tibble(dups), 
                                     dplyr::as_tibble(missing), 
                                     proportion = duplication[dupIND]) %>% 
            dplyr::mutate(silly = silly) %>% 
            dplyr::select(silly, tidyr::everything())
          
        }
        if(!sum(dupIND)){
          
          report <- tibble::tibble(silly = silly, ID1 = NA, ID2 = NA, Missing1 = NA, Missing2 = NA, proportion = NA) 
          
        }
        
      }
      
      DupDist <- tibble::tibble(silly = silly, IDs = names(duplication), duplication) %>% 
        tidyr::separate(IDs, into = c("ID1", "ID2"), sep = "[:punct:]")
      
      list(report = report, DupDist = DupDist)  
      
    }  # End multicore loop
    
    parallel::stopCluster(cl)
    
    report <- lapply(1:length(sillyvec), function(i){
      
      dupcheck[[i]][["report"]]
      
    }) %>% dplyr::bind_rows()
    
    DupDist <- lapply(1:length(sillyvec), function(i){
      
      dupcheck[[i]][["DupDist"]]
      
    }) %>% dplyr::bind_rows()
    
    output <- list(report = report, DupDist = DupDist)
    
    print(Sys.time() - start.time)
    
    return(output)
    
  }  # end if quantile
  
}