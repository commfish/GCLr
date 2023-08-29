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
#' This function uses [hierfstat::hierfstat()] to calculate locus statistics for Ho, Fis, and Fst.
#' 
#' @seealso 
#' [hierfstat::hierfstat()]
#' [hierfstat::varcomp()]
#' [hierfstat::basic.stats()]
#' [GCLr::create_hierfstat_data()]
#' 
#' @examples
#'   
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#'   
#' GCLr::locus_stats(sillyvec = sillyvec, loci = GCLr::ex_LocusControl$locusnames[-c(10, 12, 13, 32, 33)], ncores = parallel::detectCores(), LocusCtl = GCLr::ex_LocusControl )
#' 
#' @export
locus_stats <- function(sillyvec, loci, ncores = 4, LocusCtl = LocusControl) {
    
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  start.time <- Sys.time()
  
  ploidy <- LocusCtl$ploidy[loci]
  
  dat <- GCLr::create_hierfstat_data(sillyvec = sillyvec, region = NULL, pop = seq_along(sillyvec), loci = loci, ncores = ncores, LocusCtl = LocusCtl)

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