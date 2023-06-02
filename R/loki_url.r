
#' Build URL to connect to LOKI
#' 
#' This function builds the URL to connect to LOKI and is called by other GCL functions that connect to LOKI via [DBI::dbConnect()].
#'
#' @note this function is used in the background
#' 
#' @return The URL string used by [DBI::dbConnect()] to connect to the Oracle cloud database (LOKI)
#' 
#' @aliases LOKI_URL.GCL
#' 
#' @export
#' 
loki_url <- function(test.db = FALSE){


  if(test.db == TRUE){
    
    host <- "soaora6-scan.us1.ocm.s7134325.oraclecloudatcustomer.com"
    
    port <- "1521"
    
    svc <- "DFGGCLT.us1.ocm.s7134325.oraclecloudatcustomer.com"
    
    protocol <-"jdbc:oracle:thin:"
    
  }else{
    
    host <- "soaora7-scan.us1.ocm.s7134325.oraclecloudatcustomer.com"
    
    port <- "1521"
    
    svc <- "dfgcfresp.us1.ocm.s7134325.oraclecloudatcustomer.com"
    
    protocol <-"jdbc:oracle:thin:"
    
  }
  
  url <- paste0(protocol, "@", host, ":", port, "/", svc)
  
  return(url)
  
}

#' @rdname loki_url
#' @export
LOKI_URL.GCL <- loki_url  