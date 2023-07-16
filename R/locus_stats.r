#' @title Calculate Basic Locus Statistics
#' 
#' @description
#' This function calculates basic locus statistics (Ho, Fis, Fst) for each locus across collections (silly code).
#'
#' @param sillyvec A character vector of silly codes without the ".gcl" extension.
#' @param loci A character vector of locus names.
#' @param ncores A numeric value for the number of cores to use in a \pkg{foreach} `%dopar%` loop (default = 4). 
#' The number of cores cannot exceeds the number on your device ([parallel::detectCores()]).
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()] (default = LocusControl).
#'
#' @returns A tibble with 1 row per locus in `loci` + "Overall" and the following 4 columns:
#'     \itemize{
#'       \item \code{locus}: locus name
#'       \item \code{Ho}: observed heterozygosity calculated by [hierfstat::basic.stats()]
#'       \item \code{Fis}: Fis calculated by [hierfstat::varcomp()]
#'       \item \code{Fst}: Fst calculated by [hierfstat::varcomp()]
#'       }       
#'       
#' @details
#' This function uses the [hierfstat::hierfstat()] to calculate locus statistics for Ho, Fis, and Fst.
#' Previous versions of this function used [GCLr::gcl2fstat()] to write out an FSTAT (*.dat) file, but this is no longer the case.
#' Use [GCLr::gcl2fstat()] if you want to create an FSTAT (*.dat) file.
#' 
#' @seealso 
#' [hierfstat::hierfstat()]
#' [hierfstat::varcomp()]
#' [hierfstat::basic.stats()]
#' [GCLr::create_hierfstat_data()]
#' 
#' @examples
#' \dontrun{
#' GCLr::create_locuscontrol(markersuite = "Sockeye2011_96SNPs",
#'                           username = .username,
#'                           password = .password)
#' sillyvec = c("SKARLSE11", "SKARLSE99L", "SREDCRY11", "SDOGSC08")
#' GCLr::loki2r(sillyvec = sillyvec,
#'              username = .username,
#'              password = .password)
#' (my_locus_stats <- GCLr::locus_stats(sillyvec = sillyvec, loci = LocusControl$locusnames, ncores = 4))
#' }
#'  
#' @export
locus_stats <- function(sillyvec, loci, ncores = 4, LocusCtl = LocusControl, ...){

  unused_args <- list(...)
  
  if (!length(unused_args) == 0){ 
    
    if(names(unused_args) %in% c("fstatdir", "dir")){
      
      warning("This function no longer writes out an FSTAT .dat file; therefore; the 'fstatdir' and 'dir' arguemnts are no longer used. 
            Use gcl2fstat to produce an FSTAT file if needed.")
      
      }
    
    }
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  start.time <- Sys.time()
  
  ploidy <- LocusCtl$ploidy[loci]
  
  dat <- GCLr::create_hierfstat_data(sillyvec = sillyvec, region = NULL, pop = seq_along(sillyvec), loci = loci, ncores = ncores)

  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores) #Start cluster  
  
  #Get variance components
  #Start parallel loop
  
  `%dopar%` <- foreach::`%dopar%`
  
  MyVC <- foreach::foreach(locus = loci, .packages = "hierfstat") %dopar% {
    
    diploid <- ploidy[locus]==2
    
    if(diploid){
      
      VC <- hierfstat::varcomp(dat[, c("pop", locus)], diploid = diploid)$overall
      
      table <- tibble::tibble(locus = !!locus, P = VC[1], I = VC[2], G = VC[3])
      
    }
    
    if(!diploid){
      
      VC <- hierfstat::varcomp(dat[, c("pop", locus)], diploid = diploid)$overall
      
      table <- tibble::tibble(locus = !!locus, P = VC[1], I = NA, G = VC[2])
      
    } 
    
    table
    
  } %>%  dplyr::bind_rows()
  
  parallel::stopCluster(cl)# Stop cluster
  
  #Summarize variance components
  MyTable0 <- MyVC %>% 
    dplyr::rowwise() %>%  
    dplyr::mutate(total = sum(c(P,I,G), na.rm = TRUE)) %>% 
    dplyr::mutate(P = P/total, I = I/total) %>% 
    dplyr::select(locus, P, I, total) %>% 
    dplyr::ungroup()
  
  Overall <- MyVC %>% 
    dplyr::summarize(P = sum(P, na.rm = TRUE), I = sum(I, na.rm = TRUE), total = sum(MyTable0$total, na.rm = TRUE)) %>% 
    dplyr::mutate(locus = "Overall", I = I/total, P = P/total)
  
  MyTable <- dplyr::bind_rows(MyTable0, Overall) 
  
  #Heterozygosities
  Hovec <- apply(hierfstat::basic.stats(dat[, c("pop", loci[ploidy==2])])$Ho, 1, mean, na.rm = TRUE)
  
  Ho <- tibble::tibble(locus = c(loci[ploidy==2], "Overall"), Ho = c(Hovec, mean(Hovec, na.rm = TRUE)))
  
  #Join varcomp summary with Ho  
  output <- MyTable %>% 
    dplyr::left_join(Ho, by = "locus") %>% 
    dplyr::mutate(Fis = I, Fst = P) %>% 
    dplyr::select(locus, Ho, Fis, Fst)
  
  print(Sys.time() - start.time)

  return(output)
  
}