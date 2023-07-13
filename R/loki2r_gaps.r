#' Create GCL Objects (GAPS)
#'
#' This function connects to LOKI and creates a ".gcl" object for each silly code in `sillyvec`. A "*.gcl" object is a Gene Conservation Laboratory genotypes object with associated sample attributes.
#'
#' @param sillyvec Vector of silly codes to pull from LOKI.
#' @param username Username for LOKI connection.
#' @param password Password for LOKI connection.
#' 
#' @return The runtime of the function in seconds.
#' 
#' @details This function requires the R package "RJDBC" for connecting to LOKI.
#' It connects to LOKI using the provided \code{username} and \code{password}, and retrieves genotypes and sample attributes for each silly code in \code{sillyvec}.
#' The genotypes and attributes are stored in separate "*.gcl" objects named after each silly code.
#'
#' @examples
#' \dontrun{
#' sillyvec <- c("KQUART06", "KQUART08", "KQUART09")
#' loki2r_gaps(sillyvec = sillyvec, username = "awbarclay", password)
#' }
#' 
#' @export
loki2r_gaps <- function(sillyvec, username, password){

  start.time <- Sys.time() 
  
  options(java.parameters = "-Xmx10g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  markersuite <- "GAPS_Chinook_uSATs"
  
  lociqry <- paste0("SELECT * FROM AKFINADM.V_LOCUS_MARKERSUITE WHERE SUITE_NAME = '", markersuite, "'")
  
  locidata <- RJDBC::dbGetQuery(con, lociqry)
  
  locusnames <- as.character(locidata$LOCUS_NAME)
  
  Publishedlocusnames <- as.character(locidata$PUBLISHED_NAME)
  
  nloci <- length(locusnames)
  
  ploidy <- sapply(locusnames, function(locus) {
    qry <- paste0("SELECT PLOIDY, LOCUS_NAME FROM AKFINADM.LOCUS_LOOKUP WHERE LOCUS_NAME IN ", "'", locus, "'")
    2^(as.character(RJDBC::dbGetQuery(con, qry)$PLOIDY[1]) == "D")
  } )
  
  ConTable <- RJDBC::dbGetQuery(con, "SELECT * FROM AKFINADM.GAPS_ALLELE_CONVERSION")
  
  ConTable$VALUE_ADFG <- as.numeric(ConTable$VALUE_ADFG)  # this is crucial for allele sorting so that "99" is at the front, not the back of the allele list
  
  ConTable$VALUE_CTC <- as.numeric(ConTable$VALUE_CTC)  # this is crucial for allele sorting so that "99" is at the front, not the back of the allele list
  
  alleles <- sapply(locusnames, function(loc) {as.character(sort(unique(as.vector(ConTable[ConTable[, "LOCUS_NAME"] == loc, "VALUE_CTC"]))))} )
  
  nalleles <- sapply(alleles, function(allele) {length(allele)} )
  
  LocusCtl <- list(MarkerSuite = paste0(markersuite, "_Converted"), 
       locusnames = locusnames, 
       Publishedlocusnames = Publishedlocusnames, 
       alleles = alleles, 
       nalleles = nalleles, 
       ploidy = ploidy)
  
  assign("LocusControl", LocusCtl,
         pos = 1)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  loci <- LocusCtl$locusnames
  
  nloci <- length(loci)
  
  ploidy <- LocusCtl$ploidy
  
  alleles <- LocusCtl$alleles
  
  nalleles <- LocusCtl$nalleles
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(silly in sillyvec){
    
    collectionIDqry = paste0("SELECT * FROM AKFINADM.GEN_COLLECTIONS WHERE SILLY_CODE = ", "'", silly, "'")
    
    collectionID = RJDBC::dbGetQuery(con, collectionIDqry)$COLLECTION_ID
    
    gnoqry = paste0("SELECT * FROM AKFINADM.V_GAPS_GENOTYPES WHERE SUITE_NAME = '", markersuite, "' AND FK_COLLECTION_ID IN (", collectionID, ") ORDER BY LOCUS, FK_COLLECTION_ID,FK_FISH_ID")
    
    my.gno = RJDBC::dbGetQuery(con, gnoqry)
    
    if(nrow(my.gno) == 0) {
      
      message(paste0(silly, ".gcl has no data, skipping ", match(silly, sillyvec), " of ", length(sillyvec), " failed."))
      
      next()
      
      }
    
    fishID = as.character(unique(my.gno$FK_FISH_ID))
    
    nind = length(fishID)
    
    my.gno$LOCUS <- as.factor(my.gno$LOCUS)
    
    sillyloci = levels(my.gno$LOCUS)[match(loci, levels(my.gno$LOCUS))[!is.na(match(loci, levels(my.gno$LOCUS)))]]
    
    lociIND = is.na(match(loci, sillyloci))
    
    if(sum(lociIND)) {stop(paste0("No available data for loci '", loci[lociIND], "' for silly '", silly, "', hoser!!!"))}
    
    scoredIND = table(my.gno$FK_FISH_ID, my.gno$LOCUS)[fishID, loci, drop = FALSE]
    
    includeIND = apply(scoredIND, 1, sum) == nloci
    
    scores = abind::abind(tapply(X = my.gno$ALLELE1_CONV, INDEX = data.frame(my.gno$FK_FISH_ID, my.gno$LOCUS),
                          FUN = function(allele) {as.character(allele)} )[fishID, loci, drop = FALSE],
                   tapply(X = my.gno$ALLELE2_CONV, INDEX = data.frame(my.gno$FK_FISH_ID, my.gno$LOCUS),
                          FUN = function(allele) {as.character(allele)} )[fishID, loci, drop = FALSE], along = 3)
    
    dimnames(scores)[[3]] = c("Dose_1", "Dose_2")
    
    scores[scores == 0] = NA
    
    gntps = array(NA, c(nind, nloci, max(nalleles)), dimnames = list(fishID, loci, paste0("Allele ", 1:max(nalleles))))
    
    for(ind in fishID) {
      for(locus in loci) {
        for(al in 1:nalleles[locus]) {
          gntps[ind, locus, al] = sum(scores[ind, locus, 1:ploidy[locus]] == alleles[[locus]][al])
        }  # al
        if(sum(gntps[ind, locus, 1:nalleles[locus]]) != ploidy[locus] | sum(is.na(gntps[ind, locus, 1:nalleles[locus]])) > 0) {
          gntps[ind, locus ,1:nalleles[locus]] = rep(NA, nalleles[locus])
        }  # if
      }  # locus
    }  # ind
    
    attrbtqry = paste0("SELECT * FROM AKFINADM.GEN_SAMPLED_FISH_TISSUE WHERE FK_COLLECTION_ID=", collectionID)
    
    attributes = RJDBC::dbGetQuery(con, attrbtqry)
    
    if(nrow(attributes) == 0) {
      attributes = data.frame(array(NA, c(nind, ncol(attributes)), dimnames = list(fishID, dimnames(attributes)[[2]])))
      attributes$FK_FISH_ID = as.numeric(fishID)
      attributes$FK_COLLECTION_ID = rep(collectionID, nind)
    }

    plateID = sapply(fishID, function(fshID){
      qry <- paste0("SELECT FK_PLATE_ID FROM AKFINADM.GEN_DNA_WELL WHERE SILLY_CODE=", paste0("'", silly, "'"), " AND FISH_NO=", fshID)
      IDs = unlist(RJDBC::dbGetQuery(con, qry))
      if(length(IDs) == 0) {IDs=NA}
      nIDs = length(IDs)
      if(nIDs == 1) {return(IDs)}
      if(nIDs > 1) {return(paste(IDs[1:nIDs], collapse = "/"))} 
    } )
    
    IND = match(fishID, attributes$FK_FISH_ID)
    
    attributes = attributes[IND, ]
    
    dimnames(attributes)[[1]] = fishID
    
    SillySource = paste(silly, fishID, sep = "_")
    
    attributes = cbind(attributes, plateID, SillySource)
    
    assign(paste0(silly, ".gcl"), list(counts = gntps, scores = scores, n = nind, attributes = attributes), pos = 1)
    
    message(paste0(silly, ".gcl created ", match(silly, sillyvec), " of ", length(sillyvec), " completed."))
  }  # silly
  
  discon <- RJDBC::dbDisconnect(con)
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime) 
  
  return(fulltime)
}
 