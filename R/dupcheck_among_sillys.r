#' Check for Duplicate Individuals Between ".gcl" Objects
#'
#' This function checks for duplicate individuals between silly code ".gcl" objects. It is mainly used for quality control purposes.
#'
#' @param KeySillys A character vector of silly codes without the ".gcl" extension.
#' @param KeySillyIDs A named list of character vectors containing FK_FISH_IDs for each silly code in `KeySillys` to check against all silly codes in `BetweenSillys`. If NULL, all FK_FISH_IDs for each silly code in `KeySillys` are checked against all silly codes in `BetweenSillys`.
#' @param BetweenSillys A character vector of silly codes without the ".gcl" extension.
#' @param loci A character vector of locus names. If set to NULL, all loci in the ".gcl" objects will be used.
#' @param minnonmissing The proportion of loci that a pair must share non-missing genotypes in order to be reported.
#' @param minproportion The proportion of shared non-missing loci that must be shared between the individuals to be reported as a matching pair.
<<<<<<< HEAD
#' @param ncores The number of cores to use in a `foreach::%dopar%` loop. If the number of cores exceeds the number on your device, then `ncores` defaults to `parallel::detectCores()`.
=======
#' @param ncores The number of cores to use in a [foreach::%dopar%] loop. If the number of cores exceeds the number on your device, then `ncores` defaults to [parallel::detectCores()].
>>>>>>> 4e8ba2320f2cc09da272f26eff0e5386c70fce79
#' @param plot.results Logical value indicating whether to produce histograms of duplicate rates for each silly code in `KeySillyIDs`.
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)  
#'
#' @return A tibble that includes the duplicate individuals that exceeded `minproportion` or had the maximum duplicate rate, as well as duplicate individuals that came from the same original collection (i.e., project duplicates). The tibble includes the following variables:
#'   \itemize{
#'     \item Keysillyvial: The KeySilly_ID.
#'     \item Betweensillyvial: The BetweenSilly_ID.
#'     \item Keymissing: The number of loci without genotypes for each KeySilly_ID.
#'     \item Betweenmissing: The number of loci without genotypes for each BetweenSilly_ID.
#'     \item DuplicateRate: The proportion of duplicate genotypes between each KeySilly_ID and BetweenSilly_ID.
#'   }
#' The function also prints a histogram of duplicate rates for each silly code in `KeySillyIDs`. Each plot has a vertical red line indicating the `minproportion`, and the highest duplicate rate bar is labeled with the `BetweenSillyIDs` with that duplicate rate. The title of each plot indicates the `KeySillyID` that was checked for duplicates.
#'
#' @examples
#' password <- "************"
#' create_locuscontrol(markersuite = "CookInletChinook2013_43SNPs", username = "awbarclay", password = password)
#' loki2r(sillyvec = c("KKILL05", "KKILL06", "KFUNN05", "KFUNN06"), username = "awbarclay", password = password)
#' pool_collections(collections = c("KKILL05", "KFUNN05"), IDs = list(KKILL05 = KKILL05.gcl$FK_FISH_ID, KFUNN05 = c(9, 30)), newname = "KKILL05")
#' pool_collections(collections = c("KFUNN06", "KKILL06"), IDs = list(KFUNN06 = KFUNN06.gcl$FK_FISH_ID, KKILL06 = c(101, 176)), newname = "KFUNN06")
#' KKILL05.gcl <- KKILL05.gcl %>% mutate(SillySource = case_when(FK_FISH_ID %in% c(69, 70) ~ paste(SILLY_CODE, FK_FISH_ID, sep = "_"), TRUE ~ SillySource))
#' KFUNN06.gcl <- KFUNN06.gcl %>% mutate(SillySource = case_when(FK_FISH_ID %in% c(184, 185) ~ paste(SILLY_CODE, FK_FISH_ID, sep = "_"), TRUE ~ SillySource))
#' KKILL06qc.gcl <- KKILL06.gcl %>% mutate(SILLY_CODE = "KKILL06qc") %>% mutate(SillySource = paste(SILLY_CODE, FK_FISH_ID, sep ="_"))
#' KFUNN05qc.gcl <- KFUNN05.gcl %>% mutate(SILLY_CODE = "KFUNN05qc") %>% mutate(SillySource = paste(SILLY_CODE, FK_FISH_ID, sep ="_"))
#' KeySillys <- c("KKILL06qc", "KFUNN05qc")
#' KeySillyIDs <- list(KKILL06qc = c(101, 176), KFUNN05qc = c(9, 30))
#' BetweenSillys <- c("KKILL05", "KKILL06", "KFUNN05", "KFUNN06")
#' loci <- LocusControl$locusnames
#' results <- GCLr::dupcheck_among_sillys(KeySillys = c("KKILL06qc", "KFUNN05qc"), KeySillyIDs = list(KKILL06qc = c(101, 176), KFUNN05qc = c(9, 30)), BetweenSillys = c("KKILL05", "KKILL06", "KFUNN05", "KFUNN06"), loci = LocusControl$locusnames, minnonmissing = 0.6, minproportion = 0.9, ncores = 4)
#'
#' @export
dupcheck_among_sillys <- function(KeySillys, KeySillyIDs = NULL, BetweenSillys, loci, minnonmissing = 0.6, minproportion = 0.9, ncores = 4, plot.results = TRUE, LocusCtl = LocusControl){

  if(!all(loci %in% LocusCtl$locusnames)){
    
    stop(paste0("The following `loci` were not found in `LocusControl`:\n", paste(setdiff(loci, LocusCtl$locusnames), collapse = "\n")))
    
  }
  
  if(sum(KeySillys %in% BetweenSillys) > 1){
    
    stop("One or more of the KeySillys are also included in BetweenSillys. This is not the intended use of this function. Use dupcheck_within_silly() to find duplicates within a silly.")
    
    }
    
  
  if(ncores > parallel::detectCores()) {
    
    ncores <- parallel::detectCores()
    
  }
  
  ploidy <- LocusCtl$ploidy[loci]
  
  scores_cols <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
    as.vector()  # This keeps the scores columns in the correct order when there are loci with similar names.
  
  # my.key
  if(!is.null(KeySillyIDs)){
    
    my.key <- lapply(KeySillys, function(silly){
      
      get(paste0(silly, ".gcl")) %>% 
        filter(FK_FISH_ID %in% KeySillyIDs[[silly]])
      
    }) %>% bind_rows()
      
  } else{
    
    my.key <- lapply(KeySillys, function(silly){
      
      get(paste0(silly, ".gcl"))
      
    }) %>% bind_rows()
    
  }
  
  # my.between - this is a list for parallel loop
  my.between <- lapply(BetweenSillys, function(silly){
    
    get(paste0(silly, ".gcl"))
    
  }) %>% set_names(BetweenSillys)
  
  # Added this if statement code for haploid markers, rubias::close_matching_samples() was counting them as missing loci because of the NAs in the allele2 column. 
  # This can be removed Eric Anderson fixes the function. 
  if(any(ploidy == 1)) {
    
    haploci <- names(ploidy[ploidy == 1])
    
    for(locus in haploci){
      
      my.key[[paste0(locus, ".1")]] <- my.key[[locus]]
      
      for(silly in BetweenSillys){
        
         my.between[[silly]][[paste0(locus, ".1")]] <- my.between[[silly]][[locus]]
        
      }
      
    }
  } 
    
  ## Loop through between sillys
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  `%dopar%` <- foreach::`%dopar%`
  
  dupcheck0 <- foreach::foreach(silly = BetweenSillys, .export = c("loci"), .packages = c("tidyverse","rubias")) %dopar% {
    
    new.gcl <- my.between[[silly]] %>%
      dplyr::bind_rows(my.key) %>% 
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
    
    dups <- rubias::close_matching_samples(D = new.gcl, gen_start_col = 5, min_frac_non_miss = minnonmissing, min_frac_matching = 0) %>% 
      dplyr::filter(collection_1 == silly, collection_2 %in% KeySillys)
    
    dups %>% 
      dplyr::filter(collection_1 == silly, collection_2 %in% KeySillys) %>% 
      dplyr::rename(Betweensillyvial = indiv_1, Keysillyvial = indiv_2, Betweensilly = collection_1, Keysilly = collection_2) %>% 
      dplyr::mutate(DuplicateRate = num_match/num_non_miss) %>% 
      dplyr::select(Keysillyvial, Betweensillyvial, DuplicateRate)
    
  } %>% dplyr::bind_rows()  # End multicore loop
  
  parallel::stopCluster(cl) 
  
  # Plots
  
  if(plot.results == TRUE){ #Added this as an option to not produce plots
    
    sapply(dupcheck0$Keysillyvial %>% unique(), function(key){
      
      dc <- dupcheck0 %>% 
        tidyr::separate(Betweensillyvial, into = c("Betweensilly", NA), sep = "_", remove = FALSE) %>% 
        dplyr::filter(Keysillyvial == key)
      
      max_dup <- dc %>% 
        dplyr::filter(DuplicateRate == max(DuplicateRate))
      
      plot <- dc %>% 
        ggplot2::ggplot(aes(x = DuplicateRate, fill = Betweensilly)) + # I added "fill = BetweenSilly" so the histogram will show which between sillys make up the distribution - this may be removed if others don't like it.
        ggplot2::geom_histogram(bins = 100) +
        ggplot2::geom_vline(xintercept = minproportion, color = "red", size = 1.25)+
        ggplot2::xlab("Duplicate Rate")+
        ggplot2::ylab("Frequency") +
        ggplot2::xlim(0, 1.02) +
        ggplot2::ggtitle(label = paste0("KeySillyID: ", key))+
        ggplot2::geom_text(aes(x = max_dup$DuplicateRate %>% unique, y = length(max_dup$DuplicateRate)), label = paste0(max_dup$Betweensillyvial, collapse = "_"), angle = 90, hjust = -.05)
      
      suppressWarnings(print(plot))
      
    }) # End plots
    
  }
    
  max_dups <- dupcheck0 %>% 
    dplyr::filter(DuplicateRate == max(DuplicateRate))
  
  threshold_dups <- dupcheck0 %>% 
    dplyr::filter(DuplicateRate >= minproportion)
   
  project_dups <- dupcheck0 %>% 
    dplyr::mutate(Keysillyvial_2 = gsub(pattern = "qc_", replacement = "_", x = Keysillyvial)) %>% 
    dplyr::filter(Keysillyvial_2 == Betweensillyvial) %>% 
    dplyr::select(Keysillyvial, Betweensillyvial, DuplicateRate)
  
  suppressMessages(dupcheck <- dplyr::bind_rows(max_dups, threshold_dups, project_dups) %>% unique)
  
  scores1 <- paste0(loci, ".1")
  
  Keymissing <- my.key %>% 
    dplyr::select(Keysillyvial = SillySource, all_of(scores1)) %>% 
    dplyr::mutate(across(all_of(scores1), is.na)) %>% 
    dplyr::group_by(Keysillyvial) %>% 
    dplyr::mutate(Keymissing = sum(!!!syms(scores1))) %>% 
    dplyr::select(Keysillyvial, Keymissing)
  
  Betweenmissing <- my.between %>% 
    dplyr::bind_rows() %>% 
    dplyr::filter(SillySource %in% dupcheck0$Betweensillyvial) %>% 
    dplyr::select(Betweensillyvial = SillySource, all_of(scores1)) %>% 
    dplyr::mutate(across(all_of(scores1), is.na)) %>%  
    dplyr::group_by(Betweensillyvial) %>% 
    dplyr::mutate(Betweenmissing = sum(!!!syms(scores1))) %>% 
    dplyr::select(Betweensillyvial, Betweenmissing)
  
  duplicate_summary <- dupcheck %>% 
    dplyr::left_join(Keymissing, by = "Keysillyvial") %>% 
    dplyr::left_join(Betweenmissing, by = "Betweensillyvial") %>% 
    dplyr::select(Keysillyvial, Betweensillyvial, Keymissing, Betweenmissing, DuplicateRate)
  
  return(duplicate_summary)
  
}