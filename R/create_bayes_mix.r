#' @title Create BAYES mixture files
#'
#' @description This function creates a BAYES mixture (.mix) file for each silly in mixvec.
#'
#' @param mixvec Character vector of ".gcl" objects you want to produce mixture files for.
#' @param loci Character vector of the loci you wish to include.
#' @param dir Character vector of where to save the ".bse" file.
#' @param ncores Number of cores to use for parallel processing.
#'
#' @returns Returns the fortran format of the mixture file, which is needed for [GCLr::create_bayes_ctl()].
#' It also saves each mixture in `mixvec` as a .mix file.
#' 
#' @note If you want to analyze more than one silly as a mixture, use [GCLr::pool_collections()] to combine them into a new silly.gcl
#'
#' @examples
#' \dontrun{
#' load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
#' loki2r(sillyvec = c("KSUSC18FW", "KSUSCN18", "KSUSC19FW", "KSUSCN19"), username = "awbarclay", password = password)
#' combine_loci(sillyvec = c("KSUSC18FW", "KSUSCN18", "KSUSC19FW", "KSUSCN19"), markerset = c("Ots_MHC1", "Ots_MHC2"))
#' pool_collections(collections = c("KSUSC18FW", "KSUSCN18"), newname = "Susitna2018")
#' pool_collections(collections = c("KSUSC19FW", "KSUSCN19"), newname = "Susitna2019")
#' loci <- c(loci82, "Ots_MHC1.Ots_MHC2")
#' mix_fortran <- create_bayes_mix(mixvec = c("Susitna2018", "Susitna2019"), loci = loci, dir = getwd(), ncores = 8)
#' }
#' 
#' @export
create_bayes_mix <- function(mixvec, loci, dir, ncores = 4) {

    if (sum(is.na(match(loci, LocusControl$locusnames)))) {
    
    stop(paste("'", loci[is.na(match(loci, LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  start_time <- Sys.time()
  
  my.gcl <- sapply(mixvec, function(silly){
    
    get(paste(silly, ".gcl", sep = ""), pos = 1)
    
  }, simplify = FALSE)
   
  for (silly in mixvec) { #Start silly loop
    
    new.gcl <- my.gcl[[silly]]
  
    scores_names <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
      as.vector() 
    
    scores0 <- new.gcl %>% 
      dplyr::select(tidyselect::all_of(scores_names)) %>% 
      as.data.frame(stringsAsFactors = FALSE)
    
    id <- as.numeric(dimnames(scores0)[[1]])
    
    cl <- parallel::makePSOCKcluster(ncores)
    
    doParallel::registerDoParallel(cl, cores = ncores)  
    
    #Had to suppress an anoying message which only printed to console when running the foreach loop.
    suppressMessages(
      
      mix_scores <-  foreach::foreach(loc = loci, .packages = c("tidyverse", "tidyselect", "janitor"), .export = c("LocusControl"), .combine = dplyr::bind_cols) %dopar% { #Start multicore loop: loci
      
      variables <- c(loc, paste(loc, 1, sep = "."))
        
      my.alleles <- LocusControl$alleles[[loc]]
        
      comb_alleles0 <- scores0 %>%
        dplyr::select(tidyselect::all_of(variables)) %>% 
        dplyr::mutate(dplyr::across(dplyr::everything(), ~factor(.x, levels = my.alleles$call))) %>% 
        dplyr::mutate(dplyr::across(dplyr::everything(), ~as.numeric(.x))) %>% 
        tibble::as_tibble() %>% 
        dplyr::mutate(id = id) %>% 
        tidyr::pivot_longer(-id) %>% 
        dplyr::mutate(value = factor(value, levels = my.alleles$allele)) %>% 
        janitor::tabyl(id, value, show_missing_levels = TRUE, show_na = TRUE)
      
      comb_alleles <-  comb_alleles0 %>% 
        dplyr::select(dplyr::all_of(my.alleles$allele %>% as.character())) %>% 
        tidyr::unite(col = "comb", as.character(my.alleles$allele), sep = "") %>%
        dplyr::pull(comb) %>% 
        stringr::str_pad(width = length(my.alleles$allele) + 1, pad = " ", side = "left") %>% 
        tibble::as_tibble()
      
      comb_alleles
      
      } %>% 
      tidyr::unite(col = "mix", dplyr::everything(), sep = "") #End multicore loop: loci
    
      ) #Suppress
    
    write.table(mix_scores, file = paste0(dir, "/", silly, ".mix"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    parallel::stopCluster(cl)
      
    } #End silly loop
  
  a_rle <- sapply(loci, function(loc){LocusControl$nalleles[[loc]]}) %>% 
    rle()
  
  mix_fortran <- paste0("(", paste0(sapply(seq(length(a_rle$values)), function(locus){
    
    paste(a_rle$lengths[locus], "(1X,", a_rle$values[locus], "I1)")
    
    }), collapse = ","), ")")
  
  Sys.time() - start_time
  
  return(mix_fortran)
  
}