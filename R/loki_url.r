loki_url <- function(test.db = FALSE){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function builds the URL to connect to LOKI and is called on by other GCL 
  #  functions that connect to LOKI via DBI::dbConnect.
  # 
  #  Output:
  #  The URL string used by DBI::dbConnect to connect to the Oracle cloud database (LOKI)
  #
  #  Written by Andy Barclay 4/15/19
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
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
  
