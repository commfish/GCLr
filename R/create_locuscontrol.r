#' Create LocusControl Object
#'
#' This function connects to LOKI (GCL database) and creates a "LocusControl" object.
#'
#' @param markersuite optional; the pre-defined name in LOKI for the set of markers for which you want genotypes (e.g., markersuite="KenaiChinook2010_40SNPs"). This set must be pre-defined in LOKI. 
#'                   
#' @param locusnames optional; a character vector of locus names spelled exactly the way they are in LOKI. Use "Locus Report Brief" in OceanAK to get the correct locus names.
#' 
#' @param species optional; argument to find markersuite by species. (i.e. sockeye, chinook, chum, coho; default = NULL)
#' 
#' @param username your state user name
#' 
#' @param password your password used to access LOKI - contact Eric Lardizabal if you don't have a password for LOKI
#'
#' @return This function assigns a tibble with the following columns: 
#'   \itemize{
#'     \item \code{Markersuite} (character): markersuite name.
#'     \item \code{locusnames} (character): the locus name used by the GCL.
#'     \item \code{Publishedlocusnames} (character): the locus name used in publications.
#'     \item \code{nalleles} (integer): the number of alleles for each locus.
#'     \item \code{ploidy} (integer): the ploidy of the locus (2 = diploid, 1 = haploid).
#'     \item \code{alleles} (named list): named list of tibbles containing the alleles for each locus.
#'   }
#'   The tibble will be named "LocusControl" and assigned to your current workspace.
#'   
#' @details
#' There are three ways to create a LocusControl object using this function: 1) supply a markersuite name, 2) supply a list of locusnames (this will become a user-defined markersuite), or 3) supply the species name and select from a list of available markersuites in the console.
#' Only one of these three ways can be used at a time, so make sure the arguments you aren't using are NULL (see examples)
#'   
#' @examples
#' \dontrun{
#' create_locuscontrol(markersuite = "Sockeye2011_96SNPs", locusnames = NULL, species = NULL, username = "awbarclay", password = password)
#' create_locuscontrol(markersuite = NULL, locusnames = NULL, species = NULL username = "awbarclay", password = password)
#' create_locuscontrol(markersuite = NULL, locusnames = c("One_ACBP-79","One_agt-132","One_aldB-152","One_apoe-83" ), species = NULL, username = "awbarclay", password = password)
#'}
#'
#' @export
create_locuscontrol <- function(markersuite = NULL, locusnames = NULL, species = NULL, username, password){
  
  if(exists("LocusControl", where = 1)) {
    
    stop("LocusControl already exists")
    
  }
  
  if(is.null(markersuite) & is.null(locusnames) & is.null(species)){
    
    stop("The user must supply one of the following arguments: markersuite, locusnames, species")
    
  }
  
  if(!is.null(markersuite) & !is.null(locusnames) & !is.null(species)){
    
    stop("The user must supply only one of the following arguements: markersuite, locusnames, species")
    
  }
  
  if(!is.null(markersuite) & !is.null(locusnames) & is.null(species)){
    
    stop("The user must supply only one of the following arguements: markersuite, locusnames, species")
    
  }
  
  if(!is.null(markersuite) & is.null(locusnames) & !is.null(species)){
    
    stop("The user must supply only one of the following arguements: markersuite, locusnames, species")
    
  }
  
  if(is.null(markersuite) & !is.null(locusnames) & !is.null(species)){
    
    stop("The user must supply only one of the following arguements: markersuite, locusnames, species")
    
  }
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")

  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
   
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  if(!is.null(species)){
    
    lociqry <- "SELECT * FROM AKFINADM.V_LOCUSQRY"  # Query locus information of markers in markersuite.
    
    mkrsuites <- RJDBC::dbGetQuery(con, lociqry)
    
    options <- mkrsuites %>% 
      tibble::as_tibble() %>% 
      dplyr::filter(grepl(pattern = species, x = SUITE_NAME, ignore.case = TRUE)) %>% 
      dplyr::pull(SUITE_NAME) %>% 
      unique()
    
    markersuite <- select.list(options, "Select a markersuite:", multiple = FALSE)
    
  }
  
  # Query by 'markersuite', else query by 'locusnames'
 
  if(is.null(locusnames) & !is.null(markersuite)) {
    
    lociqry <- paste("SELECT * FROM AKFINADM.V_LOCUSQRY WHERE SUITE_NAME = '", markersuite, "'", sep = "")  # Query locus information of markers in markersuite.
    
    locidata <- RJDBC::dbGetQuery(con, lociqry)
    
    locusnames <- sort(locidata$LOCUS_NAME)
    
  } else {
    
    markersuite <- "User defined from locusnames"
    
    lociqry <- paste("SELECT * FROM AKFINADM.V_LOCUSQRY WHERE LOCUS_NAME IN (",paste0("'", locusnames, "'", collapse = ","), ")", sep = "")  # Query locus information of markers in locusnames.
    
    locidata <- RJDBC::dbGetQuery(con, lociqry)
    
  }
  
  # Warn user if some 'locusnames' not found
  if(!all(locusnames %in% locidata$LOCUS_NAME)) {
    
    miss_loci <- locusnames[!locusnames %in% locidata$LOCUS_NAME]
    
    nmiss_loci <- length(miss_loci)
    
    message(paste(nmiss_loci, "out of", length(locusnames), "locusnames not found in LOKI!!!"))
    
    sapply(miss_loci, function(locus) {
      
      message(locus)
      
    })
    
    locusnames <- locusnames[locusnames %in% locidata$LOCUS_NAME]
    
  }
  
  Publishedlocusnames <- purrr::set_names(locidata$PUBLISHED_NAME, locidata$LOCUS_NAME)[locusnames]
  
  nloci <- length(locusnames)
  
  ploidy <- purrr::set_names(locidata$PLOIDY, locidata$LOCUS_NAME)[locusnames]
  
  alleleqry <- paste("SELECT * FROM AKFINADM.V_ALLELEQRY WHERE LOCUS IN (", paste0("'", locusnames, "'", collapse = ","), ")", sep = "")  # Get possible alleles from allele lookup table.
  
  alleles0 <- RJDBC::dbGetQuery(con, alleleqry)
  
  discon <- RJDBC::dbDisconnect(con) #Disconnect from database
  
  alleles <- sapply(locusnames, function(locus){
    
    my.loc = dplyr::filter(alleles0, LOCUS==locus) %>% 
      dplyr::select(-VALUE_CTC)
    
    nalleles = dim(my.loc)[1] 
    
    tibble::tibble(allele = seq(nalleles), call = my.loc$VALUE %>% sort())
    
  }, simplify = FALSE)
  
  nalleles <- sapply(alleles, dim)[1,] 
  
  assign("LocusControl", tibble::tibble(MarkerSuite = markersuite, locusnames = locusnames, Publishedlocusnames = Publishedlocusnames, nalleles = nalleles, ploidy = ploidy, alleles = alleles), pos = 1, envir = .GlobalEnv) #Assign elements to LocusControl list.	
  
  ans <- paste0("LocusControl created for markersuite ='", markersuite, "'")
  
  return(ans)
  
}