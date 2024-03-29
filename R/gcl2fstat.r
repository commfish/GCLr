#' Convert ".gcl" Objects to FSTAT
#'
#' Write out an FSTAT input file from ".gcl" objects.
#'
#' @param sillyvec A character vector of silly codes to include in the FSTAT file.
#' @param loci A character vector of locus names to include in the FSTAT file.
#' @param path The full file path to write out the FSTAT file including ".dat" extension; with "\\" or "/" separator between folders.
#' @param ncores The number of cores for multithreading using [doParallel()] and [foreach()]. Default is 4. 
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)  
#' 
#' @details 
#' This function requires a `LocusControl` object. Run [GCLr::create_locuscontrol()] prior to this function. 
#' Use this function when you want to run an analysis in `FSTAT`. 
#' 
#' @return Writes out an FSTAT file to the specified path.
#' 
#' @examples
#' 
#' sillys <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' loci <- GCLr::ex_baseline[,-c(1:5)] %>%
#'   names() %>%
#'   gsub(pattern = "*\\.1", x = ., replacement = "") %>%
#'   unique()
#'   
#' GCLr::gcl2fstat(sillyvec = sillys, loci = loci, path = path.expand("~/FSTATfile.dat"), ncores = 4, LocusCtl = GCLr::ex_LocusControl)
#' 
#' @seealso [GCLr::create_hierfstat_data()]
#' 
#' @export
gcl2fstat <- function(sillyvec, loci, path, ncores = 4, LocusCtl = LocusControl){

  start_time <- Sys.time()
  
  if(sum(is.na(match(loci, LocusCtl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocusCtl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  nsillys <- length(sillyvec)
  
  nloci <- length(loci)
  
  ploidy <- LocusCtl$ploidy[loci]
  
  alleles <- LocusCtl$alleles[loci] %>% 
    dplyr::bind_rows(.id = "locus")
      
  maxchar <- max(nchar(alleles$allele))+1
    
  nalleles <- LocusCtl %>% 
    dplyr::filter(locusnames%in%loci) %>% 
    dplyr::pull(nalleles) %>% 
    purrr::set_names(loci)
  
  my.gcl <- lapply(sillyvec, function(silly){
    
    get(paste(silly, ".gcl", sep = ""), pos = 1)
    
    }) %>% 
    purrr::set_names(sillyvec)
  
  # Multicore loop
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  `%dopar%` <- foreach::`%dopar%`
  
  scores_all <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse")) %dopar% {
    
    new.gcl <- my.gcl[[silly]]
    
    pop_no <- match(silly, sillyvec)
    
    scores_names <- sapply(loci, function(locus) {
      
      c(locus, paste0(locus, ".1"))
      
      }) %>% 
      as.vector() 
    
    scores <- new.gcl %>% 
      dplyr::select(tidyselect::all_of(scores_names)) %>% 
      as.data.frame(stringsAsFactors = FALSE)
    
    lapply(loci, function(loc){
      
      variables <- c(loc, paste(loc, 1, sep = "."))
      
      my.alleles <- alleles %>% 
        dplyr::filter(locus==loc)
      
      scores %>%
        dplyr::select(tidyselect::all_of(variables)) %>% 
        dplyr::mutate(dplyr::across(dplyr::everything(), .fns = ~factor(., levels = my.alleles$call))) %>% 
        dplyr::mutate(dplyr::across(dplyr::everything(), .fns = as.numeric)) %>% 
        dplyr::mutate(dplyr::across(dplyr::everything(), .fns = as.character)) %>%
        tidyr::replace_na(replace = list("0", "0") %>%
                            set_names(variables)) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), .fns = ~stringr::str_pad(., width = maxchar, pad = "0", side = "left"))) %>% 
        tidyr::unite(col = !!rlang::as_name(loc), tidyselect::all_of(variables), sep = "")
      
    }) %>% 
      dplyr::bind_cols() %>% 
      purrr::set_names(loci) %>% 
      dplyr::mutate(Pop = pop_no) %>% 
      tidyr::unite("comb_scores", c(Pop, tidyselect::all_of(loci)), sep = " ") 
  
  } # end multicore loop
  
  parallel::stopCluster(cl)
  
  fstat <- paste(nsillys, nloci, max(nalleles), maxchar, sep = " ")
  
  fstat <- rbind(fstat, cbind(loci))
  
  fstat <- rbind(fstat,  scores_all %>% unlist() %>% cbind())
  
  write.table(x = fstat, file = path, row.names = FALSE, col.names = FALSE, quote = FALSE) 
  
  message("\nfstat file created\n", sep = '')
  
  Sys.time()-start_time
  
}