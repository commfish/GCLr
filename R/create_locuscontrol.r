create_locuscontrol <- function(markersuite = NULL, locusnames = NULL, username, password){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function connects to LOKI and creates a "LocusControl" object  
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   markersuite - "markersuite" is the pre-defined name in LOKI for the set of markers for which you want genotypes (e.g. markersuite="KenaiChinook2010_40SNPs"). 
  #                  This set must be pre-defined in LOKI (see Eric Lardizabal). 
  #
  #   locusnames - a character vector of locus names spelled exactly the way they are in LOKI.
  #                Use Locus Report Brief in Ocean AK to get the correct locus names.
  #
  #   username - your state user name
  #
  #   password - your password used to access LOKI - see Eric Lardizabal if you don't have a passord for LOKI
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function assigns a tibble with the following columns: Markersuite <character>: markersuite name, 
  #                                                              locusnames <character>: the locus name used by the GCL, 
  #                                                              Publishedlocusnames<character>: the locus name used in publications,
  #                                                              nalleles <int>: the number of alleles for each locus,
  #                                                              ploidy <int>: the ploidy of the locus (2 = diploid, 1 = haploid),
  #                                                              alleles <named list>: named list of tibbles containing the alleles for each locus
  #   The tibble will be named "LocusControl"
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   password = "************"
  # 
  #   create_locuscontrol(markersuite = "UCI_Chinook_GTSeq_557SNPs", locusnames = NULL, username = "awbarclay", password = password)
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  if(exists("LocusControl", where = 1)){
    
    stop("LocusControl already exists")
    
  }
  
  
  # This copies the "odbc8.jar" file to the R folder on your computer if it doesn't exist there. This file contains the java odbc drivers needed for RJDBC
  
  if(!file.exists(path.expand("~/R"))){ 
    
    dir <- path.expand("~/R")
    
    dir.create(dir)
    
    bool <- file.copy(from ="V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
    
  } else {
    
    if(!file.exists(path.expand("~/R/ojdbc8.jar"))){
      
      bool <- file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
      
    }
    
  }
  
  options(java.parameters = "-Xmx10g")
  
  if(file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
    
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc8.jar", " ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
  } else {
    
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = path.expand("~/R/ojdbc8.jar"), " ")
    
  } 
  
  url <- loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  # Query by 'markersuite', else query by 'locusnames'
  
  if(is.null(markersuite) & is.null(locusnames)) {stop("Need to provide either 'locusnames' or 'markersuite'")}
  
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