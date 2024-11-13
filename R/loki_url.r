#' Build URL to connect to LOKI
#' 
#' This function builds the URL to connect to LOKI and is called by other GCL functions that connect to LOKI via [DBI::dbConnect()].
#'
#' @note this function is used in the background
#' 
#' @return The URL string used by [DBI::dbConnect()] to connect to the Oracle cloud database (LOKI)
#' 
#' @keywords internal
loki_url <- function(test.db = FALSE){

  if(test.db == TRUE){
    
    host <- "soaocinp-db6hq-scan.exa.sjcprod.oraclevcn.com"
    
    port <- "1521"
    
    svc <- "DFGGCLT.exa.sjcprod.oraclevcn.com"
    
    protocol <-"jdbc:oracle:thin:"
    
  }else{
    
    host <- "soaocip-ovpnc-scan.exa.sjcprod.oraclevcn.com"
    
    port <- "1521"
    
    svc <- "DFGGCLP.exa.sjcprod.oraclevcn.com"
    
    protocol <-"jdbc:oracle:thin:"
    
  }
  
  url <- paste0(protocol, "@", host, ":", port, "/", svc)
  
  return(url)
  
}