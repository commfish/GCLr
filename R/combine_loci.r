#' Combine Loci
#'
#' This function combines a set of markers into a single marker. 
#'
#' @param sillyvec A character vector of silly codes.
#' @param markerset A vector of the set of loci you wish to combine.
#' @param update Logical switch. If TRUE, the "LocusControl" object is updated and all "*.gcl" objects in "sillyvec" will be updated with the new marker. If FALSE, the "LocusControl" object is not updated and a temporary object called "*.temp.gcl" with the updated data is created.
#' @param delim Specifies the separator between combined loci, either a period (.) which is the default or an underscore (_) so locus names will work in SPAM.
#'
#' @examples
#' sillyvec <- c("KQUART06","KQUART08","KQUART09")
#' markerset <- c("Ots_vatf-251", "Ots_ZR-575")
#' combine_loci(sillyvec, markerset, update = TRUE, delim = c(".", "_")[1])
#' combine_loci(sillyvec, markerset, update = FALSE, delim = c(".", "_")[1])
#'
#' @details
#' This function requires a LocusControl object. Run [GCLr::create_locuscontrol()] prior to this function. This function requires dplyr version 1.0.0 or higher.
#'
#' @aliases CombineLoci.GCL
#'
#' @export

combine_loci <- function(sillyvec, markerset, update = TRUE, delim = c(".", "_")[1]){

  if(!all(markerset %in% LocusControl$locusnames)){
    
    stop(paste0("'", setdiff(markerset, LocusControl$locusnames), "' from argument 'markerset' not found in 'LocusControl' object!!!"))
    
  }
  
  nmarkers <- length(markerset)  
  
  myploidy <- LocusControl$ploidy[markerset]
  
  if(sum(myploidy == myploidy[1]) != nmarkers){
    
    stop("'markerset' has different ploidies!!!")
    
  }  
  
  MarkerSuite <- LocusControl$MarkerSuite %>% 
    unique()
  
  locusnames <- LocusControl$locusnames
  
  newmarkername <- paste(markerset, collapse = delim)
  
  existnewmarker <- newmarkername %in% locusnames
  
  loci <- unique(c(locusnames, newmarkername))
  
  nloci <- length(loci)
  
  Publishedlocusnames <- LocusControl$Publishedlocusnames
  
  Publishedlocusnames <- c(Publishedlocusnames, purrr::set_names(NA, newmarkername))
  
  newalleles <- get_phenotypes(markerset)
  
  maxchar <- max(nchar(newalleles))
  
  alleles <- LocusControl$alleles[locusnames]
  
  nalleles <- LocusControl$nalleles[locusnames]
  
  ploidy <- LocusControl$ploidy[locusnames]
  
  if(!existnewmarker){
    
    alleles[[length(locusnames) + 1]] <- tibble::tibble(allele = seq_along(newalleles), call = newalleles)
    
    nalleles <- c(nalleles, length(newalleles)) 
    
    ploidy <- c(ploidy, 1) 
    
  }  
  
  names(alleles) <- names(nalleles) <- names(ploidy) <- loci
  
  for(silly in sillyvec){
    
    my.gcl <- get(paste0(silly, ".gcl"), pos = 1)
    
    if(newmarkername %in% names(my.gcl)){
      
      warning(paste0("'", newmarkername, "'"," already created in silly '", silly,"'!!!"))
      next()
      
    }
    
    newmarkername_1 <- paste0(newmarkername, ".1")  # Allele 2 name
    
    # Combine haploid
    if(unique(myploidy) == 1){ 
      
      new.gcl <- my.gcl %>% 
        tidyr::unite(col = {{newmarkername}}, tidyselect::all_of(markerset), sep = '', remove = FALSE, na.rm = TRUE) %>%  # Had to add the {{}} around the col object for this to work. 
        mutate(!!rlang::sym(newmarkername) := dplyr::case_when(nchar(!!rlang::sym(newmarkername)) < maxchar ~ NA_character_,  # Added this mutate case_when to replace any genotypes with less than maxchar with NA's. Maybe there is a better way?
                                                               TRUE ~ !!rlang::sym(newmarkername)), 
               !!rlang::sym(newmarkername_1) := NA_character_) %>% 
        dplyr::relocate(!!rlang::sym(newmarkername), .after = tidyselect::last_col()) %>% 
        dplyr::relocate(!!rlang::sym(newmarkername_1), .after = tidyselect::last_col())  # Had to relocate alleles separately or they get reordered
      
    }
    
    # Combine diploid
    if(unique(myploidy)==2){ 
      
      sel_var <- lapply(markerset, function(mkr){
        c(mkr, paste0(mkr, ".1"))
        }) %>% unlist() #Setting up order of variable to unite so the markers are in the same order as markerset
      
      new.gcl <- my.gcl %>% 
        tidyr::unite(col = {{newmarkername}}, tidyselect::all_of(sel_var), sep = '', remove = FALSE, na.rm = TRUE) %>%  # Had to add the {{}} around the col object for this to work. 
        dplyr::mutate(!!rlang::sym(newmarkername) := dplyr::case_when(nchar(!!rlang::sym(newmarkername)) < maxchar ~ NA_character_,  # Added this mutate case_when to replace any genotypes with less than maxchar with NA's. Maybe there is a better way?
                                                               TRUE ~ !!rlang::sym(newmarkername)), 
               !!rlang::sym(newmarkername_1) := NA_character_) %>% 
        dplyr::relocate(!!rlang::sym(newmarkername), .after = tidyselect::last_col()) %>% 
        dplyr::relocate(!!rlang::sym(newmarkername_1), .after = tidyselect::last_col())  # Had to relocate alleles separately or they get reordered
      
    }
    
    if(update){
      
      assign(paste0(silly, ".gcl"), new.gcl, pos = 1)
      
    }
    
    if(!update){
      
      assign(paste0(silly, ".temp.gcl"), new.gcl, pos = 1)
      
    } 
    
  }
  
  if(update){
    
    assign(
      x = "LocusControl",
      value = tibble::tibble(
        MarkerSuite = MarkerSuite,
        locusnames = loci,
        Publishedlocusnames = Publishedlocusnames,
        nalleles = nalleles,
        ploidy = ploidy,
        alleles = alleles
      ),
      pos = 1
    )
    
  }   
  
}
#' @rdname combine_loci
#' @export
CombineLoci.GCL <- combine_loci  