#' Create rubias Mixture
#'
#' @description This function creates the mixture dataframe needed for \pkg{rubias}.
#'
#' @param sillyvec A character vector of the mixture silly codes, each silly code is treated as its own mixture.
#' @param loci A character vector of the loci to include.
#' @param path A character vector of where to save each mixture as a ".csv" (see details).
#' 
#' @details 
#' This function saves each mixture as its own ".csv" file for posterity. To read in these ".csv" files, use \code{readr::read_csv(file = file, col_types = cols(.default = "c")} to make sure all columns are character vectors (if homozygous for T, it will become a logical vector).
#' 
#' @return A dataframe in \pkg{rubias} mixture format and writes out a ".csv" file of each mixture.
#' 
#' @seealso [GCLr::ex_mixtures]
#' 
#' @examples
#' \dontrun{
#'  
#'  mixvec <- c("S1_2018", "S2_2018", "S3_2018", "S4_2018", "S5_2018")
#'  
#'  loci <- LocusControl$locusnames
#'  
#'  GCLr::create_rubias_mix(sillyvec = mixvec, loci = loci, path = "rubias/mixture")
#'  
#' }
#' 
#' @export
create_rubias_mix <- function(sillyvec, loci, path = "rubias/mixture") {

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