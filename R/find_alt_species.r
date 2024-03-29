#' Identify Alternate Species
#'
#' This function identifies the likelihood of incorrect species based on a few identifiable markers. Specifically, it works with chum and sockeye.
#'
#' @param sillyvec A character vector of silly codes (e.g., c("SYEHR07", "SWENA98", "STATS05")).
#' @param species The species you are working with; options = sockeye, chum; default = chum.
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)  
#' 
#' @details 
#' This function requires a `LocusControl` object and ."gcl" objects. Run [GCLr::create_locuscontrol[] and [GCLr::loki2r] prior to this function.
#' The Alternate and Failed marker lists were text files in a random folder on a random network drive; marker lists were hard-coded into function so it doesn't relying on hidden files. In the future, if applicable, we can simply add more markers here.
#' 
#' @return A plot of Failed vs Alternate and a tibble containing:
#'   \itemize{
#'     \item silly_fish: The fish ID.
#'     \item prop_alter: The proportion of alternate markers.
#'     \item prop_fail: The proportion of failed markers.
#'     \item num_alter: The number of non-missing alternate markers used for analyses.
#'   }
#'   
#' @examples
#' \dontrun{
#' 
#' wrong_spp <- GCLr::find_alt_species(sillyvec = sillys, species = "sockeye")
#' 
#' }
#' 
#' @export
find_alt_species <- function(sillyvec, species = "chum", LocusCtl = LocusControl){

  # First, load the alternate & failed markers for a given species, or stop if incorrect species is provided. 
  if (species == "chum") {
    
    AlternateGenotypes <- dplyr::tribble(~AlternateMarker,	~AlternateGenotype,
                                         "Oke_U1018-50", "TT",
                                         "Oke_U2032-74",	"GG",
                                         "Oke_Cr386",	"GH")

    FailedMarkers <- dplyr::tribble(~FailedMarker,
                                 "Oke_AhR1-78",
                                 "Oke_CKS1-94",
                                 "Oke_e2ig5-50",
                                 "Oke_U1002-262",
                                 "Oke_U1025-135",
                                 "Oke_u200-385",
                                 "Oke_U2025-86",
                                 "Oke_U502-241")
  
    } else if (species == "sockeye") {
    
    AlternateGenotypes <- dplyr::tribble(~AlternateMarker, ~AlternateGenotype,
                                         "One_ctgf-301",	"TT",
                                         "One_KCT1-453",	"GG",
                                         "One_taf12-248",	"TT",
                                         "One_U1013-108",	"GG",
                                         "One_U1214-107",	"AA")
    
    FailedMarkers <- dplyr::tribble(~FailedMarker,
                                    "One_Hsp47",
                                    "One_STC-410",
                                    "One_U1010-81",
                                    "One_MHC2_190",
                                    "One_cin-177",
                                    "One_vatf-214")
    } else {
      
      stop( "This can only be run on sockeye or chum.")
  }

  loci <- LocusCtl$locusnames # List of loci from LocusControl; Relies on LocusControl in environment!


  # Subset alternate markers by only those found in LocusControl
  AlternateGenotypes <- AlternateGenotypes %>% 
    dplyr::filter(AlternateMarker %in% loci) 
  
  # How many alternate markers are present in Locus Control?
  AlternateCount <- AlternateGenotypes %>% 
    dplyr::tally() %>% 
    dplyr::pull(n)
  
  if ( AlternateCount == 0 ) { # If no alternate markers found in LocusControl,
    
    stop("None of the alternate markers found in LocusControl.") # Halt - no matching markers
    
  } else { # Otherwise, subset based on which alternate markers are present
    
    message(
      
      paste0("Performing analyses on ", AlternateCount, " out of ", length(AlternateGenotypes$AlternateMarker), " alternative markers.\n") # Pasting in values to display number of markers
     
      ) # Warning msg displaying number of alternate markers found in LocusControl
    
    }
 
  # Subset failed markers by only those found in LocusControl
  MyFailedMarkers <- FailedMarkers %>% 
    dplyr::filter(FailedMarker %in% loci) %>% 
    dplyr::pull(FailedMarker)
  
  # How many failed markers are present in Locus Control?
  FailedCount <- length(MyFailedMarkers) # Total number of failed markers found in LocusControl
  
  if ( FailedCount == 0 ) { # If no failed markers found in LocusControl,
    
    stop("None of the Failed markers found in LocusControl.") # Halt - no matching markers
    
  } else { # Otherwise, subset based on which failed markers are present
    
    message(
      
      paste0("Performing analyses on ", FailedCount, " out of ", length(FailedMarkers$FailedMarker), " failed markers.\n") # Pasting in values to display number of markers
    
      ) # Warning msg displaying number of failed markers found in LocusControl
    
  }
  
  alternate_vars <- c(AlternateGenotypes$AlternateMarker, paste0(AlternateGenotypes$AlternateMarker, ".1")) %>%
    sort() # Vector of variable names to select alternate and failed loci explicitly.
  
  failed_vars <- c(MyFailedMarkers, paste0(MyFailedMarkers, ".1")) %>%
    sort() # Vector of variable names to select alternate and failed loci explicitly. 
  
  # Create tibble of all fish in sillyvec, containing just alternate and failed loci
  gclobjectsAll <- sapply(sillyvec, function(silly){
    
    get(paste0(silly, ".gcl")) %>% # Grab gcl objects
      dplyr::select(SILLY_CODE, FK_FISH_ID, dplyr::all_of(c(alternate_vars, failed_vars))) %>% # Select fish ids and subset by markers of intrest
      tidyr::unite(silly_fish, SILLY_CODE, FK_FISH_ID) # Make fish_ID col for convenience
    
  }, simplify = FALSE) %>% # Keep nested
    dplyr::bind_rows() # Turn nested list into tibble
 
 # Calculate failure
  Failure <- gclobjectsAll %>%
    dplyr::select(c(silly_fish, dplyr::all_of(MyFailedMarkers))) %>%  # Filter by just failed markers
    dplyr::mutate(failure = 
                    rowSums(is.na(dplyr::select(., -silly_fish))) / FailedCount) %>% # Calculate fail
    dplyr::select(silly_fish, failure) # Grab important columns

 # Create tibble of alternates
  gclobjectsAllAlternate <- lapply(AlternateGenotypes$AlternateMarker, function(locus){
    
    locus_vars <-  c(locus, paste0(locus, ".1")) # List of loci, both versions
    
    gclobjectsAll %>% 
      tidyr::unite(!!dplyr::sym(locus), !!!dplyr::syms(locus_vars), sep = "") %>% # Merge locus and locus.1 calls
      dplyr::select(!!dplyr::sym(locus)) # Select just the locus
    
  }) %>% 
    dplyr::bind_cols() %>% # Adding locus from apply into tibble cols
    dplyr::mutate(silly_fish = gclobjectsAll$silly_fish) %>% # Pull in ID column
    dplyr::select(silly_fish, dplyr::everything()) %>% # Sort, IDs first followed by loci
    dplyr::mutate(dplyr::across(dplyr::all_of(AlternateGenotypes$AlternateMarker), 
                                ~gsub(x = .x, pattern = "NANA", replacement = NA) # Replace any NANA which were created by joining locus+locus.1 containing NAs
                  ))
  
 # Now calculations to find alternate spp.
  Alternate <- gclobjectsAllAlternate %>%
    tidyr::pivot_longer(-silly_fish, names_to = "locus", values_to = "call") %>% # Make long format
    dplyr::left_join(AlternateGenotypes, by = c("locus" = "AlternateMarker")) %>% # Join the list of alternates
    dplyr::mutate(GenoMatch = 
                    dplyr::case_when(
                      is.na(call)               ~ NA_real_, # If genotype is NA then assign NA
                      call == AlternateGenotype ~ 1, # Does genotype equal alternate genotype?
                      TRUE                      ~ 0)) %>% # Otherwise, genotype is not the alternate
    dplyr::group_by(silly_fish) %>% 
    dplyr::summarize(alternate = sum(GenoMatch, na.rm = TRUE)/AlternateCount, # Calculate prop. alternates
              non_missing_alt = sum(!is.na(call)), # Calculate number of non-NA alternate markers (i.e., # out of total available for ea. fish)
              .groups = "drop_last") 

   Results <- Failure %>% dplyr::left_join(Alternate, by = "silly_fish") # Tibble with failure, alternate, and non_missing_alternate calculations
  
   plotly::ggplotly(
     Results %>% 
       dplyr::filter(!is.na(alternate)) %>%
       ggplot2::ggplot(ggplot2::aes(x = failure, y = alternate)) +
       ggplot2::geom_count(alpha = 0.5) +
       ggplot2::labs(title = "Proportion of failure vs alternate indicator markers",
            x = "Failure", 
            y = "Alternate", 
            size = 12)
   ) %>% print() # Trick to print plotly plots AND return tibble of results! 
   
 return(Results)
  
}

