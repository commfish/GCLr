#' @title Calculate Basic Locus Statistics
#' 
#' @description
#' This function calculates basic locus statistics for each locus across collections (silly code).
#' 
#' @param data a hierfstat data object (default = NULL)
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
#'       \item \code{Hs}: gene diversity calculated by [hierfstat::basic.stats()]
#'       \item \code{Ar}: allelic richness calculated by [hierfstat::allelic.richness()]
#'       \item \code{Fis}: Fis calculated by [hierfstat::varcomp()]
#'       \item \code{Fst}: Fst calculated by [hierfstat::varcomp()]
#'       }       
#'       
#' @details
#' This function uses \pkg{hierfstat} to calculate locus statistics for Ho, Hs, Ar, Fis, and Fst. If a hierfstat data object is supplied (recommended), sillyvec and loci objects
#' are not needed, and the function runs a lot faster.  If a hierfstat data object is not supplied (i.e., data = NULL) and sillyvec and loci objects are supplied, the
#' function will create a temporary hierfstat data object to calculate the locus statistics. The hierfstat data object produced by this function is not assigned to your workspace, so it's best 
#' to create a data object using [GCLr::create_hierfstat_data()] if you want to use the object for something else.
#' 
#' @seealso 
#' [hierfstat::varcomp()]
#' [hierfstat::basic.stats()]
#' [hierfstat::allelic.richness()]
#' [GCLr::create_hierfstat_data()]
#' 
#' @examples
#'   
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' loci <- ((get(paste0(sillyvec[1], ".gcl")) %>% names())[-c(1:19)])[seq(1, 185, by = 2)]
#' 
#' dat <- GCLr::create_hierfstat_data(sillyvec, region = NULL, pop = seq_along(sillyvec), loci = loci, ncores = parallel::detectCores(), LocusCtl = GCLr::ex_LocusControl)
#' 
#' GCLr::locus_stats(data = dat, sillyvec = NULL, loci = NULL, ncores = parallel::detectCores(), LocusCtl = GCLr::ex_LocusControl)#Supplying data object
#'
#' GCLr::locus_stats(data = NULL, sillyvec = sillyvec, loci = loci, ncores = parallel::detectCores(), LocusCtl = GCLr::ex_LocusControl)#Supplying sillyvec and loci
#' 
#' @export
locus_stats <- function(data = NULL, sillyvec = NULL, loci = NULL, ncores = parallel::detectCores(), LocusCtl = LocusControl) {
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  if(is.null(data) & is.null(sillyvec) & is.null(loci)){
    
    stop("Please supply a hierfstat data object or sillyvec and loci objects.")
    
    }
    
  start.time <- Sys.time()
  
  if(!is.null(data)){
    
    dat <- data
    
    loci <- data %>% 
      dplyr::select(-any_of(c("region", "pop", "spop"))) %>% 
      names()
    
  } else(
    
    if(is.null(sillyvec)|is.null(loci)){
      
      stop("sillyvec and loci must be supplied when data = NULL.")
      
      }else(dat <- GCLr::create_hierfstat_data(sillyvec = sillyvec, region = NULL, pop = seq_along(sillyvec), loci = loci, ncores = ncores, LocusCtl = LocusCtl))
    
    )
 
  ploidy <- LocusCtl$ploidy[loci]
  
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
  
  #Observed Heterozygosity
  Hovec <- apply(hierfstat::basic.stats(dat[, c("pop", loci[ploidy==2])])$Ho, 1, mean, na.rm = TRUE) 
  
  Ho <- tibble::tibble(locus = c(loci[ploidy==2], "Overall"), Ho = c(Hovec, mean(Hovec, na.rm = TRUE)))
  
  #Gene diversity
  
  if(any(ploidy==2)){ # Diploid
    
    Hsvec_dip <- hierfstat::basic.stats(dat[, c("pop", loci[ploidy == 2])])$perloc[loci[ploidy == 2], ]$Hs %>% purrr::set_names(loci[ploidy==2])
  
  }else{Hsvec_dip <- NULL}
  
  if(any(ploidy==1)){ # Haploid
    
    Hsvec_hap <-  hierfstat::basic.stats(dat[, c("pop", loci[ploidy == 1])], diploid = FALSE)$perloc[loci[ploidy == 1], ]$Hs %>% purrr::set_names(loci[ploidy==1])
    
  }else{Hsvec_hap <- NULL}
  
  Hsvec <- c(Hsvec_dip, Hsvec_hap)[loci]
  
  Hs <- tibble::tibble(locus = c(loci, "Overall"), Hs = c(Hsvec, mean(Hsvec, na.rm = TRUE)))
  
  #Allelic richness
  
  if(any(ploidy==2)){ # Diploid
    
    Arvec_dip <- hierfstat::allelic.richness(dat[, c("pop", loci[ploidy == 2])])$Ar %>% rowMeans( na.rm = TRUE)
    
  }else{Arvec_dip <- NULL}
  
  if(any(ploidy==1)){ # Haploid
    
    Arvec_hap <- hierfstat::allelic.richness(dat[, c("pop", loci[ploidy == 1])])$Ar %>% rowMeans( na.rm = TRUE)
    
  }else{Arvec_hap <- NULL}
  
  Arvec <- c(Arvec_dip, Arvec_hap)[loci]
  
  Ar <- tibble::tibble(locus = c(loci, "Overall"), Ar = c(Arvec, mean(Arvec, na.rm = TRUE)))
  
  #Join varcomp summary with Ho, Hs, and Ar
  output <- MyTable %>% 
    dplyr::left_join(Ho, by = "locus") %>% 
    dplyr::left_join(Hs, by = "locus") %>% 
    dplyr::left_join(Ar, by = "locus") %>% 
    dplyr::mutate(Fis = I, Fst = P) %>% 
    dplyr::select(locus, Ho, Hs, Ar, Fis, Fst)
  
  print(Sys.time() - start.time)

  return(output)
  
}