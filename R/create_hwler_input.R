#' Create a HWLER input file
#' 
#' This function creates a HWLER input file containing mixture data only, baseline data only, or a combination mixture and baseline.
#' 
#' @param sillyvec a vector of silly codes without the ".gcl" extension. These can be mixture sillys, baseline sillys, or mixture and baseline sillys (see details)
#' 
#' @param loci character vector of locus names
#' 
#' @param samplevec a vector of numbers the same length as sillys, with baseline pops getting numbered 1:npops and mixtures getting zeros. This lets HWLER know if an individual is mixture or baseline.
#' 
#' @param file_name the path to write out the input file with .input extension (e.g. "HWLER/input/myhwler.input").
#' 
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)   
#' 
#' @returns this function returns an object called `inputfortran` containing the FORTRAN format of the input file - this object is needed for [GCLr::create_hwler_input()]
#' 
#' @details 
#' The HWLER input file can contain only mixture individuals, only baseline individuals, or a combination of mixture and baseline individuals depending on the type of analysis you are doing. 
#' 
#' @seealso See HWLER manual for additional details: [HWLER manual](system.file("HWLER", "HWLER_manual.doc", package = "GCLr"))
#' 
#' @examples
#' \dontrun{
#' 
#'  load("V:/Analysis/2_Central/Sockeye/Cook Inlet/Missing Baseline Analysis/Missing Baseline Analysis.RData")
#'  create_HWLER_mixture.GCL(sillyvec = c(PopNames69, "Western07"), loci = loci34, file_name = "HWLER/input/myHWLER.input", samplevec = c(seq(length(PopNames69)), 0))
#' 
#' }
#' 
#' @export 
create_hwler_input <- function(sillyvec, loci, file_name, samplevec, LocCtl = LocusControl){
  
  if(sum(is.na(match(loci, LocCtl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocCtl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  start_time <- Sys.time()
  
  my.gcl <- lapply(sillyvec, function(silly){
    
    get(paste(silly, ".gcl", sep = ""), pos = 1)
    
  }) %>% dplyr::bind_rows()
  
  sample_nums <- tibble::tibble(pop_num = samplevec, silly = sillyvec)
  
  scores_names <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
    as.vector() 
  
  variables <- c(loci, paste(loci, 1, sep = ".")) %>% 
    sort()
  
  my.alleles <- LocusControl$alleles %>% 
    dplyr::bind_rows(.id = "locus")
  
  max_nchar = max(nchar(samplevec))
  
  scores0 <- my.gcl %>% 
    dplyr::select(tidyselect::all_of(scores_names), silly = SILLY_CODE, ID = FK_FISH_ID) %>% 
    dplyr::left_join(sample_nums, by = "silly") %>% 
    tidyr::pivot_longer(tidyselect::all_of(variables), names_to = "locus_vars", values_to = "call") %>% 
    dplyr::mutate(locus = gsub(pattern = "\\.1", replacement = "", locus_vars )) %>% 
    dplyr::left_join(my.alleles, by = c("locus", "call")) %>% 
    dplyr::mutate(silly = factor(silly, levels = sillyvec)) %>% 
    dplyr::group_by(silly, ID, pop_num, locus) %>% 
    dplyr::summarize(count = paste0(sum(allele == 1), allele2 = sum(allele == 2)), .groups = "drop") %>% 
    dplyr::mutate(count = gsub(pattern = "NANA", x = count, replacement = "00")) %>% 
    tidyr::pivot_wider(names_from = locus, values_from = count) %>%
    tidyr::unite("silly_id", silly, ID, sep = "_")
  
  sillyID <- scores0$silly_id
  
  file <- scores0 %>% 
    dplyr::mutate(pop_num = as.character(pop_num) %>% stringr::str_pad(width = max_nchar, pad = " ", side = "left")) %>% 
    tidyr::unite(file, pop_num, tidyselect::all_of(loci), silly_id, sep = " ") %>% 
    dplyr::pull(file)
  
  a_rle <- sapply(loci, function(loc){
    
    LocusControl$nalleles[[loc]]
    
  }) %>% 
    rle()
  
  inputfortran <- paste("(","I", max_nchar,",", paste(sapply(seq(length(a_rle$values)), function(locus){
    
    paste(a_rle$lengths[locus],"(1X,", a_rle$values[locus],"I1)", sep = '')
    
  }), collapse = ","),",1X,A", max(nchar(sillyID)), ")", sep = '')
  
  write.table(file, file = file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  return(inputfortran)
  
}