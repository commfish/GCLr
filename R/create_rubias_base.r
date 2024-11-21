#' Create rubias Baseline
#'
#' This function creates the baseline dataframe needed for \pkg{rubias}.
#'
#' @param sillyvec A character vector of the silly codes for the populations to include in the baseline.
#' @param loci A character vector of the loci to include in the baseline.
#' @param group_names A character vector of reporting group names, where \code{length(group_names) == max(`groupvec`)}.
#' @param groupvec A numeric vector indicating the reporting group affiliation of each population in `sillyvec`, where \code{length(groupvec) == length(sillyvec)}.
#' @param file A character vector of where to save the baseline as ".csv" or ".fst" file (see details).
#' @param baseline_name A character vector of what to name the baseline file (see details).
#' @param out_file_type The file type for the output files, either "fst" or "csv" (default: "csv")(see details).
#' 
#' @details 
#' This function saves the baseline as a ".csv" or ".fst" file for posterity. .fst files are compressed and speed up the reading and writing process; however, they cannot be visually inspected. To read in the ".csv" files, use \code{readr::read_csv(file = file, col_types = cols(.default = "c")} to make sure all columns are character vectors (if homozygous for T, it will become a logical vector).
#' 
#' @return A dataframe in \pkg{rubias} baseline format and writes out a ".csv" or ".fst" file of the baseline.
#' 
#' @examples
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' groupvec <- GCLr::ex_baseline %>%
#'    dplyr::group_by(collection) %>%
#'    dplyr::summarize(group = unique(repunit)) %>%
#'    dplyr::mutate(group = factor(group) %>% as.numeric()) %>%
#'    dplyr::pull(group)
#' 
#' loci <- GCLr::ex_baseline[,-c(1:5)] %>%
#'    names() %>%
#'    gsub(pattern = "*\\.1", x = ., replacement = "") %>%
#'    unique()
#' 
#' group_names <- GCLr::ex_baseline$repunit %>% unique()
#' 
#' GCLr::create_rubias_base(sillyvec = sillyvec, loci = loci, group_names = group_names, groupvec = groupvec, file =  path.expand("~"), baseline_name = "my_baseline")#saves tree to documents folder
#' 
#' @export
create_rubias_base <- function(sillyvec, loci, group_names, groupvec, file = "rubias/baseline", baseline_name, out_file_type = c("fst", "csv")[2]) {

  if(!dir.exists(file)) {stop("`file` to save baseline does not exist!!!")}
  
  scores_cols <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
    as.vector() 
  
  baseline <- lapply(sillyvec, function(silly){
    
    s <- match(silly, sillyvec)
    
    get(paste0(silly, ".gcl")) %>%
      dplyr::mutate(
        sample_type = "reference",
        repunit = group_names[groupvec[s]],
        collection = silly,
        indiv = as.character(SillySource)
      ) %>%
      dplyr::select(sample_type,
                    repunit,
                    collection,
                    indiv,
                    tidyselect::all_of(scores_cols))
    
  }) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(dplyr::across(dplyr::where(is.character), ~dplyr::na_if(., "0")))  # I think this can be removed now that we convert all zeros to NAs when using loki2r
  
  if(out_file_type == "csv"){
    
    readr::write_csv(x = baseline, file = paste0(file, "/", baseline_name, "_base.csv"))
    
  }else{
    
    fst::write_fst(x = baseline, path = paste0(file, "/", baseline_name, "_base.fst"))
    
  }
  
  return(baseline)
  
}