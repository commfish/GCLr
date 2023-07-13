#' Compute global hierarchical F-statistics and variance components.
#'
#' This is a wrapper function for [hierfstat::varcomp.glob()]. This wrapper will compute variance components and hierarchical F-statistics for haploid and diploid loci and combine the results into a single output.
#'
#' @param levels A data frame containing the different levels (factors) from the outermost (e.g., region) to the innermost before the individual.
#' 
#' @param genotypes A data frame containing the genotypes in single-column numeric format (e.g., 101, 102, 202).
#' 
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()]
#'
#' @return A list containing 3 tibbles: 
#'     \itemize{
#'       \item \code{loci} - variance components by locus,  
#'       \item \code{overall} - variance components summed over all loci, 
#'       \item \code{F} - hierarchical F-statistics (see [hierfstat::varcomp.glob()] documentation).
#'     }
#' 
#' @examples
#' newbase <- GCLr::baseline %>% tidyr::separate(indiv, into = c("collection", NA), remove = FALSE)
#' 
#' sillyvec <- GCLr::base2gcl(newbase)
#' 
#' pop <- newbase %>%
#'   dplyr::group_by(collection) %>%
#'   dplyr::filter(dplyr::row_number()==1) %>%
#'   dplyr::pull(repunit) %>%
#'   factor() %>%
#'   as.numeric()
#' 
#' region <- newbase %>%
#'   dplyr::group_by(collection) %>%
#'   dplyr::filter(dplyr::row_number()==1) %>% 
#'   dplyr::mutate(region = dplyr::case_when(repunit == "KenaiOther"~1,
#'                                           TRUE~2)) %>%
#'   dplyr::pull(region)
#' 
#' loci <- GCLr::baseline[,-c(1:5)] %>%
#'   names() %>%
#'   gsub(pattern = "*\\.1", x = ., replacement = "") %>%
#'   unique()
#' 
#' fstat.dat <- GCLr::create_hierfstat_data(sillyvec = sillyvec, region = region, pop = pop, loci = loci, ncores = 4, LocusControl = GCLr::ex_LocusControl) %>% 
#'   dplyr::mutate(region = factor(region), pop = factor(pop), spop = factor(spop))
#' 
#' hierfstat_global_stats(levels = fstat.dat[ ,1:3], genotypes = fstat.dat[ ,-c(1:3)])
#'
#' @export
hierfstat_global_stats <- function(levels, genotypes, LocusCtl = LocusControl){

  loci <- names(genotypes)
  
  if(sum(is.na(match(loci, LocusCtl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocusCtl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  diploid_loci <- LocusCtl %>% 
    dplyr::filter(ploidy ==2 & locusnames %in% loci) %>% 
    dplyr::pull(locusnames)

  haploid_loci <- LocusCtl %>%
    dplyr::filter(ploidy ==1 & locusnames %in% loci) %>% 
    dplyr::pull(locusnames)

  if(length(diploid_loci) >= 1){
    
    if(length(diploid_loci) == 1){# Duplicate if only one locus, so it will work with varcomp.glob; this duplicate will be removed from the output.
      
      my.loci <- cbind(genotypes[, diploid_loci, drop = FALSE], genotypes[, diploid_loci, drop = FALSE])
      
    }else{my.loci = genotypes[, diploid_loci, drop = FALSE]}
    
    
    diploid_out <- hierfstat::varcomp.glob(levels = levels, loci = genotypes[, diploid_loci], diploid = TRUE)
      
    diploid_sum <- tibble::as_tibble(diploid_out$loc, rownames = "Locus", .name_repair = "universal") %>%
      dplyr::distinct() %>% 
      purrr::set_names(c("Locus", names(diploid_out$overall))) 
    
  }
  
  if(length(haploid_loci) >= 1){
    
    if(length(haploid_loci) == 1){# Duplicate if only one locus, so it will work with varcomp.glob; this duplicate will be removed from the output.
      
      my.loci <- cbind(genotypes[, haploid_loci, drop = FALSE], genotypes[, haploid_loci, drop = FALSE])
      
      }else{my.loci = genotypes[, haploid_loci, drop = FALSE]}
    
    haploid_out <- hierfstat::varcomp.glob(levels = levels, loci = my.loci, diploid = FALSE)
    
    haploid_sum <- tibble::as_tibble(haploid_out$loc, rownames = "Locus", .name_repair = "universal") %>%
      dplyr::distinct() %>% 
      purrr::set_names(c("Locus", names(haploid_out$overall))) 
    
  }
  
  # Combine diploid and haploid locus tables
  if(length(diploid_loci) >= 1 & length(haploid_loci) >= 1){
    
    loc_out <- dplyr::bind_rows(diploid_sum, haploid_sum)
    
  }
  
  # Only diploid locus table
  if(length(diploid_loci) >= 1 & !length(haploid_loci) >= 1){
    
    loc_out <- diploid_sum
    
  }
  
  # Only haploid locus table
  if(!length(diploid_loci) >= 1 & length(haploid_loci) >= 1){
    
    loc_out <- haploid_sum
    
  }
  
  # Sum over all loci
  overall <- loc_out %>% 
    dplyr::select(-Locus) %>% 
    apply(2, sum, na.rm = TRUE) %>% 
    t() %>% 
    tibble::as_tibble()
  
  # F-statistics matrix
  nblevels <- length(overall)
  
  fnames <- names(overall)[1:(nblevels-1)]

  f <- matrix(rep(0, (nblevels - 1)^2), ncol = (nblevels - 1))
  
  for (i in 1:(nblevels - 1)) {
    
    for (j in i:(nblevels - 1)) {
      
      f[i, j] <- sum(overall[i:j])/sum(overall[i:nblevels])
      
    }
    
  }
  
 F <- f %>%
   tibble::as_tibble(.name_repair = "universal") %>% 
   purrr::set_names(fnames)
  
output <- list(loci = loc_out, overall = overall, F = F)                      
         
return(output)                    

}