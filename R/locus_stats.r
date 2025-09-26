#' @title Calculate Basic Locus Statistics
#' 
#' @description
#' This function calculates basic locus statistics for each locus across collections (silly code).
#' 
#' @param data a hierfstat data object (default = NULL).
#' @param sillyvec A character vector of silly codes without the ".gcl" extension (default = NULL).
#' @param loci A character vector of locus names (default = NULL).
#' @param ncores ncores A numeric value for the number of cores to use in a \pkg{foreach} `%dopar%` loop (default = parallel::detectCores())
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()] (default = LocusControl).
#'
#' @returns A tibble with 1 row per locus in `loci` + "Overall" and the following 4 columns:
#'     \itemize{
#'       \item \code{locus}: locus name
#'       \item \code{Ho}: observed heterozygosity calculated by [hierfstat::basic.stats()]
#'       \item \code{Hs}: gene diversity (a.k.a. expected heterozygosity) calculated by [hierfstat::basic.stats()]
#'       \item \code{Ar}: allelic richness calculated by [hierfstat::allelic.richness()]
#'       \item \code{Fis}: Weir and Cockerham (1984) Fis calculated by [hierfstat::wc()]
#'       \item \code{Fst}: Weir and Cockerham (1984) Fst calculated by [hierfstat::wc()]
#'       }       
#'       
#' @details
#' This function uses \pkg{hierfstat} to calculate locus statistics for Ho, Hs, Ar, Fis, and Fst. If a hierfstat data object is supplied (recommended), sillyvec and loci objects
#' are not needed, and the function runs a lot faster.  If a hierfstat data object is not supplied (i.e., data = NULL) and sillyvec and loci objects are supplied, the
#' function will create a temporary hierfstat data object to calculate the locus statistics. The hierfstat data object produced by this function is not assigned to your workspace, so it's best 
#' to create a data object using [GCLr::create_hierfstat_data()] if you want to use the object for something else.
#' 
#' @seealso 
#' [hierfstat::wc()]
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
  
  ##### Weir and Cockerham (1984) Fst and Fis. Note: this has to be done for each locus separately to extract the variance components needed for calculating overall FST
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores) #Start cluster  
  
  `%dopar%` <- foreach::`%dopar%`
  
  WC_per.loc <- foreach::foreach(locus = loci, .packages = c("hierfstat", "magrittr", "tibble", "dplyr")) %dopar% {
    
    p <- ploidy[locus]==2
    
    loc.wc <- hierfstat::wc(dat[, c("pop", locus)], diploid = p)
    
    var_comps <- loc.wc$sigma %>% 
      tibble::as_tibble() %>% 
      dplyr::filter(loc == 1)
    
    a_comp <- ifelse(p == TRUE, var_comps$siga[1], sum(var_comps$siga))# Among populations variance
    
    b_comp <- ifelse(p == TRUE, var_comps$sigb[1], sum(var_comps$sigb))# Among individuals within populations variance
  
    c_comp <- ifelse(p == TRUE, var_comps$sigw[1], 0) # Within-individual variance; Haploid loci don't have this component
    
    tibble::tibble(locus = locus, Fst = loc.wc$FST, Fis = loc.wc$FIS, a = a_comp, b = b_comp, c = c_comp) %>% 
      dplyr::mutate(total_var = sum(a, b, c, na.rm = TRUE)) 
    
  } %>% dplyr::bind_rows()
  
  parallel::stopCluster(cl)# Stop cluster
  
  #Calculate overall WC Fst and Fis
  hap.loci <- loci[ploidy == 1] #FIS cannot be calculated for haploid loci, so those will be excluded from the output for that calculation
  
  Fis_overall <- WC_per.loc %>% 
    dplyr::filter(!locus %in% hap.loci) %>% 
    dplyr::summarise(Fis = sum(b)/(sum(b)+sum(c)))
    
  WC_overall <- WC_per.loc %>% 
    dplyr::summarise(locus = "Overall", Fst = sum(a)/sum(total_var)) %>% 
    dplyr::mutate(Fis = Fis_overall$Fis)
  
  WC <- dplyr::bind_rows(WC_per.loc %>% dplyr::select(locus, Fis, Fst), WC_overall)
  
  ##### Observed Heterozygosity and Expected Heterozygosity (aka Gene Diversity)
  if(any(ploidy==2)){ # Diploid
    
    dip.loci <- loci[ploidy == 2]
    
    H_dip <- hierfstat::basic.stats(dat[, c("pop", dip.loci)], diploid = TRUE)$perloc %>% 
      tibble::as_tibble(rownames = "locus") %>% 
      dplyr::select(locus, Ho, Hs)
  
  }else{H_dip <- NULL}
  
  if(any(ploidy==1)){ # Haploid
    
    hap.loci <- loci[ploidy == 1]
    
    H_hap <- hierfstat::basic.stats(dat[, c("pop", hap.loci)], diploid = FALSE)$perloc %>% 
      tibble::as_tibble(rownames = "locus") %>% 
      dplyr::select(locus, Ho, Hs)
    
    if(length(hap.loci)== 1){ # The structure of basic.stats output is different for a single haploid locus, so have treat those differently.
    
      H_hap <- H_hap[1,] %>% 
        dplyr::mutate(locus = hap.loci)
      
    }
    
  }else{H_hap <- NULL}
  
  H0 <- dplyr::bind_rows(H_dip, H_hap) 
  
  H_overall <- H0 %>% 
    dplyr::summarize(locus = "Overall", Ho = mean(Ho, na.rm = TRUE),  Hs = mean(Hs, na.rm = TRUE))
  
  H <- dplyr::bind_rows(H0, H_overall)
  
  ##### Allelic richness
  if(any(ploidy==2)){ # Diploid
    
    Arvec_dip <- hierfstat::allelic.richness(dat[, c("pop", loci[ploidy == 2])], diploid = TRUE)$Ar %>% 
      rowMeans( na.rm = TRUE)
    
  }else{Arvec_dip <- NULL}
  
  if(any(ploidy==1)){ # Haploid
    
    hap.loci <- loci[ploidy == 1]
    
    Arvec_hap <- hierfstat::allelic.richness(dat[, c("pop", hap.loci)], diploid = FALSE)$Ar %>% 
      rowMeans(na.rm = TRUE)
    
    # When there is a single haploid locus, hierfstat::allelic.richness inserts a dummy.loc as a placeholder and also replaces any dashes in the locus name with periods.
    # This if statement drops the dummy.loc and sets the correct locus name.
    if(length(hap.loci)==1){ 
      
      Arvec_hap <- purrr::set_names(Arvec_hap[1], hap.loci) 
      
    }
    
  }else{Arvec_hap <- NULL}
  
  Arvec <- c(Arvec_dip, Arvec_hap)[loci]
  
  Ar <- tibble::tibble(locus = c(loci, "Overall"), Ar = c(Arvec, mean(Arvec, na.rm = TRUE)))
  
  ##### Join varcomp summary with Ho, Hs, and Ar
  output <- WC %>% 
    dplyr::left_join(H, by = "locus") %>% 
    dplyr::left_join(Ar, by = "locus") %>% 
    dplyr::select(locus, Ho, Hs, Ar, Fis, Fst)
  
  print(Sys.time() - start.time)

  return(output)
  
}