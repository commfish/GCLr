#' Convert ".gcl" Objects to NEXUS File
#'
#' Write out a NEXUS file in GDA input format from ".gcl" objects.
#'
#' @param sillyvec A character vector of silly codes to include in the NEXUS file.
#' @param loci A character vector of locus names to include in the NEXUS file.
#' @param path The full file path to write out the NEXUS file; with "\\" or "/" separator between folders.
#' @param VialNums Logical; if TRUE (default), vial numbers will be included for each individual next to their silly code separated by an underscore (e.g., "KCRESC10_1"). If FALSE, only the SILLY code will be included for each individual (e.g., "KCRESC10").
#' @param PopNames Optional character vector the same length as `sillyvec` to give populations new names. If NULL, `PopNames` defaults to `sillyvec`.
#' @param ncores The number of cores for multithreading using [doParallel()] and [foreach()]. Default is 4. 
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)  
#' 
#' @details 
#' This function requires a `LocusControl` object. Run [GCLr::create_locuscontrol()] prior to this function.
#' 
#' @return Writes out a NEXUS file to the specified path.
#' 
#' @examples
#' sillys <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' loci <- GCLr::ex_baseline[,-c(1:5)] %>%
#'   names() %>%
#'   gsub(pattern = "*\\.1", x = ., replacement = "") %>%
#'   unique()
#' 
#' gcl2nexus(sillyvec = sillys, loci = loci, path = path.expand("~/nexusfile.nex"), VialNums = TRUE, PopNames = NULL, ncores = parallel::detectCores(), LocusCtl = GCLr::ex_LocusControl)
#' 
#' @export
gcl2nexus <- function(sillyvec, loci, path, VialNums = TRUE, PopNames = NULL, ncores = 4, LocusCtl = LocusControl){

  start_time <- Sys.time()
  
  if(sum(is.na(match(loci, LocusCtl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocusCtl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
    }
  
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  if(is.null(PopNames)) {
    
    PopNames = sillyvec
    
    }
  
  PopNames <- gsub(pattern = " ", x = PopNames, replacement = "_")
  
  ploidy <- LocusCtl %>% 
    dplyr::filter(locusnames%in%loci) %>% 
    dplyr::pull(ploidy) %>% 
    purrr::set_names(loci) %>% 
    sort(decreasing = TRUE)
  
  loci <- names(ploidy)
  
  alleles <- LocusCtl$alleles[loci] %>% 
    dplyr::bind_rows(.id = "locus")
  
  nalleles <- LocusCtl %>% 
    dplyr::filter(locusnames%in%loci) %>% 
    dplyr::pull(nalleles) %>% 
    purrr::set_names(loci)
  
  nloci <- length(loci)
  
  my.gcl <- lapply(sillyvec, function(silly){
    
    get(paste(silly, ".gcl", sep=""), pos = 1)
    
  }) %>% purrr::set_names(sillyvec)
  
  if(is.null(PopNames)){PopNames <- sillyvec}
  
  if(!length(PopNames)==length(sillyvec)){stop("PopNames is not the same length as sillyvec.")}
  
  PopNames <- PopNames %>% 
    purrr::set_names(sillyvec)
  
  file <- "#nexus"
  
  file <- rbind(file,"")
  
  file <- rbind(file, "begin gdadata; [!GDA input format]")
  
  file <- rbind(file, paste0("dimensions npops=",length(sillyvec), "  nloci=", nloci, ";"))
  
  file <- rbind(file, paste0("format missing=? separator=/;"))
  
  if(min(ploidy)==1){
    
    hapset0 <- grep("1", ploidy)
    
    if(length(hapset0)>1){hapset <- paste(range(hapset0), collapse="-")}else{hapset = hapset0}
    
    file <- rbind(file, paste0("hapset ", hapset, ";"))
    
  } else{hapset = FALSE}
  
  file <- rbind(file, "locusallelelabels")
  
  max.char <- max(nchar(seq(nloci)))
  
  locuslabels <- cbind(c(paste0("   ", format(seq(nloci-1), width = max.char), "  ", loci[seq(nloci-1)], ","),
                         paste0("   ", format(nloci, width = max.char), "  ", loci[nloci]))
                       )
  
  file <- rbind(file, locuslabels)
  
  file <- rbind(file, cbind(c(";", "matrix")))
  
  my.gcl <- lapply(sillyvec, function(silly){
    
    get(paste(silly, ".gcl", sep=""), pos = 1)
    
  }) %>% purrr::set_names(sillyvec)
  
  if(!hapset==FALSE){
    
    scores_names <- sapply(loci[-hapset0], function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
      as.vector() %>% 
      c(loci[hapset0])
    
  } else{
    
    scores_names <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
      as.vector() 
    
  }
  
  # Create scores tables
  
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  `%dopar%` <- foreach::`%dopar%`
  
  scores_all <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse", "tidyselect")) %dopar% {
    
    new.gcl <- my.gcl[[silly]]
    
    IDs <- paste(new.gcl$SILLY_CODE, new.gcl$FK_FISH_ID, sep = "_") 
    
    if(!VialNums){ 
      
      vials <- new.gcl$SILLY_CODE
      
    } else{vials <- IDs}
 
    scores <- new.gcl %>% 
      dplyr::select(tidyselect::all_of(scores_names)) %>% 
      as.data.frame(stringsAsFactors = FALSE)
    
    dimnames(scores)[[1]] = IDs 
    
    scores[scores==0] <- NA
    
    pop_scores <- lapply(loci, function(loc){ 
      
      if(ploidy[loc]==2){
        
        variables <- c(loc, paste(loc, 1, sep = "."))
                       
      } else{
                         
        variables <- loc
        
        }
      
      my.alleles <- alleles %>% 
        filter(locus==loc) %>% 
        pull(call)
      
      scores %>%
        dplyr::select(tidyselect::all_of(variables)) %>% 
        replace(is.na(.), "?") %>%
        tidyr::unite(col = loc, all_of(variables), sep = "/")
      
    }) %>% 
      dplyr::bind_cols() %>% 
      purrr::set_names(loci) %>% 
      dplyr::mutate(ID = paste0(vials)) %>% 
      tidyr::unite("comb_scores", tidyselect::all_of(loci), sep = "   ") %>% 
      tidyr::unite("ID_scores", c("ID","comb_scores"), sep = "  ") %>% 
      dplyr::pull(ID_scores) 
    
    if(silly == dplyr::last(sillyvec)){
      
      rbind(paste0(PopNames[silly],":"), 
            cbind(pop_scores), 
            t(cbind(";", "end;")))
      
    } else{
      
      rbind(paste0(PopNames[silly],":"), 
            cbind(pop_scores), 
            t(cbind("", ",")))
      
      }
    
  } #End multicore loop
  
  parallel::stopCluster(cl)
  
  file <- rbind(file, scores_all %>% unlist() %>% cbind())
  
  write.table(file, path, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  Sys.time()-start_time
 
}