#' @title Calculate individual heterozygosity
#' 
#' @description 
#' This function calculates individual heterozygosity across all non-`NA` diploid loci. This can be used as a proxy for detecting contaminated samples.
#' 
#' @param sillyvec A character vector of silly codes without the ".gcl" extension.
#' @param loci A character vector of locus names; default is all `LocusControl$locusnames`.
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()] (default = LocusControl).
#' 
#' @returns A tibble with 1 row per individual and the following 6 columns:
#'     \itemize{
#'       \item \code{SillySource}: SillySource (i.e., silly_fishID)
#'       \item \code{silly}: silly code
#'       \item \code{fish_id}: GCL fish ID
#'       \item \code{nloci}: the number of non-`NA` diploid loci genotyped per individual
#'       \item \code{ploci}: the proportion of non-`NA` diploid loci genotyped per individual
#'       \item \code{het}: individual heterozygosity (i.e., the proportion of non-`NA` diploid loci with heterozygous genotypes)
#'     }
#'   
#' @examples
#' \dontrun{
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' loci <- names(GCLr::ex_baseline)[-c(1:4)][c(TRUE, FALSE)]
#' 
#' GCLr::calc_ind_het(sillyvec = sillyvec, loci = loci, LocusCtl = GCLr::ex_LocusControl)
#' }
#' 
#' @export
calc_ind_het <- 
  function(sillyvec, 
           loci = LocusControl$locusnames, 
           LocusCtl = LocusControl) {
    
  start.time <- Sys.time() 
  
  if(!all(loci %in% LocusCtl$locusnames)){
    
    stop(paste0("'", setdiff(loci, LocusCtl$locusnames), "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
  }
  
  if(any(LocusCtl$ploidy[loci] != 2)) {
    
    message(paste0("Dropping non-diploid loci:\n", paste(names(which(LocusCtl$ploidy[loci] != 2)), collapse = "\n")))
    
    loci <- setdiff(loci, names(which(LocusCtl$ploidy[loci] != 2)))
    
  }
  
  all.gcl <- sapply(sillyvec, function(silly) {get(paste0(silly, ".gcl"), pos = 1)}, simplify = FALSE)
  
  # Start parallel loop
  het <- lapply(sillyvec, function(silly){
    
    my.gcl <- all.gcl[[silly]]
    
    dose1 <- my.gcl %>%
      dplyr::select(SillySource, dplyr::all_of(loci)) %>%
      tidyr::pivot_longer(
        cols = -SillySource,
        names_to = "locus",
        values_to = "dose1",
        values_drop_na = TRUE  # this drops any no-calls, which are NAs
      )
    
    dose2 <- my.gcl %>%
      dplyr::select(SillySource, dplyr::all_of(paste0(loci, ".1"))) %>%
      dplyr::rename_at(dplyr::vars(paste0(loci, ".1")), ~ loci) %>%
      tidyr::pivot_longer(
        cols = -SillySource,
        names_to = "locus",
        values_to = "dose2",
        values_drop_na = TRUE  # this drops any no-calls, which are NAs
      )
    
    if(all.equal(dose1$SillySource, dose2$SillySource) & all.equal(dose1$locus, dose2$locus)) {
      
      my.gcl_tall <- dplyr::bind_cols(dose1, dplyr::select(dose2, dose2))  # faster
      
    } else {
      
      my.gcl_tall <- dplyr::left_join(x = dose1, y = dose2, by = c("SillySource", "locus"))  # safer
      
    }
    
    my.gcl_tall %>%
      dplyr::group_by(SillySource) %>%
      dplyr::summarise(nloci = length(dose1),
                       ploci = nloci / length(loci),
                       het = sum(dose1 != dose2) / nloci)
    
  }) %>%
    dplyr::bind_rows()
  
  output <- het %>%
    tidyr::separate(
      col = SillySource,
      into = c("silly", "fish_id"),
      sep = "_",
      remove = FALSE
    ) %>% 
    dplyr::mutate(fish_id = as.numeric(fish_id)) %>% 
    dplyr::arrange(silly, fish_id)
  
  print(Sys.time() - start.time)
  
  return(output)
}