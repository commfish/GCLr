#' @title Calculate Allele Frequency for Multiple Loci and Silly Codes
#'
#' @description This function calculates the allele frequency for each locus for each collection (silly code).
#'
#' @param sillyvec A character vector of silly codes without the ".gcl" extension.
#' @param loci A character vector of locus names; default is all `LocusControl$locusnames`.
#' @param ncores A numeric value for the number of cores to use in a \pkg{foreach} `%dopar%` loop. 
#' If the number of cores exceeds the number on your device, `ncores` defaults to [parallel::detectCores()].
#'
#' @return A tibble with the following columns:
#'     \itemize{
#'       \item \code{silly}: silly code
#'       \item \code{locus}: locus name
#'       \item \code{allele_no}: the allele number from `LocusControl` (1 or 2 for SNPs)
#'       \item \code{allele}: the allele call from `LocusControl` (e.g. -, A, C, G, T)
#'       \item \code{freq}: the allele frequency (count)
#'       \item \code{proportion}: the allele frequency (proportion); these sum to 1 within a locus
#'     }
#'
#' @examples
#' 
#' load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
#' old2new_locuscontrol()
#' old2new_gcl(sillyvec67)
#' Freq <- calc_freq_pop(sillyvec = sillyvec67, loci = loci413)
#'
#' @export
calc_freq_pop <- function(sillyvec, loci = LocusControl$locusnames, ncores = 4){
  
  start.time <- Sys.time() 
  
  if(!all(loci %in% LocusCtl$locusnames)){
    
    stop(paste0("'", setdiff(loci, LocusCtl$locusnames), "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
  }
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  all.gcl <- sapply(sillyvec, function(silly){get(paste0(silly, ".gcl"), pos = 1)}, simplify = FALSE)
  
  scores_cols <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
    as.vector() 
  
  alleles <- LocusCtl$alleles[loci] %>% 
    dplyr::bind_rows(.id = "locus") %>% 
    dplyr::rename(allele_no = allele)
  
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  # Start parallel loop
  
  `%dopar%` <- foreach::`%dopar%`
  
  freqs <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse")) %dopar% {
    
    my.gcl <- all.gcl[[silly]]
    
    my.gcl %>% 
      dplyr::select(all_of(scores_cols)) %>% 
      tidyr::pivot_longer(cols = tidyselect::everything(), names_to = "locus", values_to = "allele") %>% 
      dplyr::mutate(locus = stringr::str_replace(string = locus, pattern = "\\.1$", replacement = "")) %>% 
      dplyr::count(locus, allele) %>% 
      dplyr::rename(freq = n) %>% 
      dplyr::right_join(alleles, by = c("locus" = "locus", "allele" = "call"), keep = TRUE) %>%
      tidyr::replace_na(list(freq = 0)) %>% 
      dplyr::mutate(silly = !!silly) %>% 
      dplyr::select(silly, locus = locus.y, allele_no, allele = call, freq) %>% 
      dplyr::arrange(locus, allele)
    
  } %>% 
    dplyr::bind_rows() 
  
  parallel::stopCluster(cl)  # End parallel loop
  
  output <- freqs %>% 
    dplyr::group_by(silly, locus) %>% 
    dplyr::mutate(total = sum(freq)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(proportion = freq/total) %>% 
    dplyr::select(-total)
  
  print(Sys.time() - start.time)
  
  return(output)
  
}