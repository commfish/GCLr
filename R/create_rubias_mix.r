create_rubias_mix <- function(sillyvec, loci, path = "rubias/mixture") {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates the mixture dataframe needed for `rubias`.
  # Each silly is treated as its own mixture.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - character vector of mixture sillys, each silly is treated as its own mixture
  #   loci - character vector of the loci you wish to include
  #   path - character vector of where to save each mixture as a .csv
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Returns a dataframe in `rubias` mixture format
  #   Saves each mixture as its own .csv file for posterity, names are derived from sillyvec
  #     NOTE: to read in these .csv's use `read_csv(file = file, col_types = cols(.default = "c"))` 
  #     to make sure all columns are character vectors (if homozygous for T, it will become a logical vector).
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # attach("V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Mixture/UCI_sockeye_2019_postseason/UCI_2019_sockeye_postseason_analysis.Rdata")
  # old2new_locuscontrol()
  # sapply(mixvec, function(mix){assign(paste0(mix, ".gcl"), get(paste0(mix, ".gcl")), pos = -1, envir = .GlobalEnv)})
  # mixvec <- mixvec
  # loci <- loci
  # detach()
  # old2new_gcl(mixvec)
  # UCI_sockeye_2019.rubias_mix <- create_rubias_mix(sillyvec = mixvec, loci = loci, path = "rubias/mixture")
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(!dir.exists(path)) {stop("`path` to save mixtures does not exist!!!")}
  
  scores_cols <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
    as.vector() 
  
  mixture <- lapply(sillyvec, function(silly){
    
    s <- match(silly, sillyvec)
    
    silly_mix <- get(paste0(silly, ".gcl")) %>%
      dplyr::mutate(
        sample_type = "mixture",
        repunit = NA_character_ ,
        collection = silly,
        indiv = SillySource
      ) %>%
      dplyr::select(sample_type,
                    repunit,
                    collection,
                    indiv,
                    tidyselect::all_of(scores_cols)) %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.character), ~dplyr::na_if(., "0")))  # I think this can be removed now that we convert all zeros to NAs when using loki2r
    
    readr::write_csv(x = silly_mix, file = paste0(path, "/", silly, "_mix.csv"))
    
  }) %>% 
    dplyr::bind_rows() 
  
  return(mixture)
  
}