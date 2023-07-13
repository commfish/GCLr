#' @title Create Pairwise Fst Tree
#'
#' @description This function generates a matrix of pairwise Fst values, a neighbor joining tree, bootstrap values for tree nodes, and variance components. 
#' It can utilize multicore processing with \pkg{foreach} to speed up the calculation of variance components.
#'
#' @param sillyvec A character vector of silly codes without the ".gcl" extension.
#' @param loci A character vector of locus names.
#' @param dir Directory to save the `PairwiseFstTree` list object using [base::dput()].
#' @param nboots A numeric value indicating the number of bootstrap iterations.
#' @param ncores A numeric value for the number of cores to use in a \pkg{foreach} `%dopar%` loop. 
#' If the number of cores exceeds the number on your device, `ncores` defaults to [parallel::detectCores()].
#' @param returnbootstrapFst A logical value indicating whether to return the Fst matrix for each bootstrap iteration.
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)  
#'
#' @return A list with the following 4 or 5 components:
#'     \itemize{
#'       \item \code{tree}: a neighbor joining tree list of 4 created by [ape::nj()] containing:
#'         \itemize{
#'           \item \code{edge}: numeric matrix
#'           \item \code{edge.length}: numeric vector
#'           \item \code{tip.lable}: character vector of population names, inherited from `sillyvec`
#'           \item \code{Nnode}: integer
#'         }
#'       \item \code{bootstrap}: numeric vector of node boodstrap values
#'       \item \code{PairwiseFst}: numeric matrix (`length(sillyvec)` x `length(sillyvec)`) of pairwise Fst values
#'       \item \code{vc}: list of `choose(n = length(sillyvec), k = 2)` with pairwise variance components from [hierfstat::varcomp()]
#'       \item \code{BootstrapFst}: list of `nboots` with numeric matrix of pairwise Fst values for each bootstrap iteration 
#'       (optional, only returned if `returnbootstrapFst = TRUE`)
#'     }
#' This list object is saved in `dir` using [base::dput()], and is named `paste(dir,"\\", length(sillyvec), "Pops", length(loci), "Loci_", "PairwiseFstTree.txt", sep = "")`
#'
#' @details Older versions of this function used to write out an FSTAT .dat file, but this is no longer the case. 
#' The function now calls on [GCLr::create_hierfstat_data()] to create a \pkg{hierfstat} data object internally.
#' Depending on the size of your baseline and number of bootstrap iterations, this function can take a while.
#'
#' @examples
#' 
#' source(paste0(path.expand("~/R/"), "Functions.R")) # GCL functions
#' load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
#' PWFSTtree <- create_pwfst_tree(sillyvec = sillyvec31, loci = loci82, dir = "output", nboots = 1000, ncores = 8, returnbootstrapFst = FALSE)
#'
#' @export
create_pwfst_tree <- function(sillyvec, loci, dir, nboots = 1000, ncores = 4, returnbootstrapFst = FALSE, LocusCtl = LocusControl){
  
   if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
  
  }
  
  message("\n4 Main function tasks: \n1) Create hierfstat data object\n2) Calculate variance componenets for each pair of sillys\n3) Calculate bootstrap Fst values\n4) Bootstrap tree nodes\n")
  
  if (.Platform$OS.type == "windows"){utils::flush.console()}
  
  start.time <- Sys.time() 
  
  n <- GCLr::silly_n(sillyvec) %>% 
    dplyr::pull(n) %>% 
    purrr::set_names(sillyvec)
  
  nsillys <- length(sillyvec)
  
  nloci <- length(loci)
  
  dat <- GCLr::create_hierfstat_data(sillyvec = sillyvec, region = NULL, pop = 1:length(sillyvec), loci = loci, ncores = ncores)

  pairs <- combn(sillyvec, 2)

  pairnames <- apply(pairs, 2, function(col){
    
    paste(col, collapse = ".")
    
    })
  
  dimnames(pairs)[[2]] <- pairnames
  
  ploidy <- LocusCtl$ploidy
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Calculate variance components for each pair of sillys
  message("\nCalculate variance components for each pair of sillys\n", sep = '')
  
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  

  `%dopar%` <- foreach::`%dopar%`
  
  vc <- foreach::foreach(pair = pairnames, .packages = "hierfstat") %dopar% {
    
    sillys <- pairs[, pair]
    
    sillyIND <- sapply(sillys, function(silly){match(silly, sillyvec)})
    
    mydat <- dat[dat[, 1]==sillyIND[1] | dat[, 1]==sillyIND[2], loci]
    
    myn <- n[sillys]
    
    levels <- rep(1:2, myn)
    
    mito.loci <- which(ploidy[loci] == 1)
    
    # Create matrix of variance components with order loci
    if(length(mito.loci) == 0) {
      
      # Diploid loci
      vc.pair <- t(sapply(loci, function(locus) {
        
        mydata = data.frame(levels, mydat[,locus])
        
        mydata_no_na <- dplyr::filter(.data = mydata, !is.na(mydat...locus.))
        
        if (dplyr::n_distinct(mydata_no_na$levels) == 2) {
          
          hierfstat::varcomp(data = data.matrix(mydata), diploid = TRUE)$overall
          
        } else {
          
          c(NA, NA, NA)
          
        }
        
      }))
      
    } else {
      
      # Diploid loci
      vc.pair <- rbind(t(sapply(loci[-mito.loci], function(locus){
        
        mydata = data.frame(levels,mydat[,locus])

        mydata_no_na <- dplyr::filter(.data = mydata, !is.na(mydat...locus.))
        
        if (dplyr::n_distinct(mydata_no_na$levels) == 2) {
          
          hierfstat::varcomp(data = data.matrix(mydata), diploid = TRUE)$overall
          
        } else {
          
          c(NA, NA, NA)
          
        }
                
      })),
      
      # Haploid loci
      t(sapply(loci[mito.loci], function(locus) {
        
        mydata <- data.frame(levels,mydat[,locus])
        
        mydata_no_na <- dplyr::filter(.data = mydata, !is.na(mydat...locus.))
        
        if (dplyr::n_distinct(mydata_no_na$levels) == 2) {
          
          c(hierfstat::varcomp(data = data.matrix(mydata), diploid = FALSE)$overall,0)
          
        } else {
          
          c(NA, NA, NA)
          
        }
        
      } ))
      
      )[loci, ]
      
    }
    
    # Replace NAs with 0s
    vc.pair[is.na(vc.pair)] <- 0
    vc.pair
  
  }

  names(vc) <- pairnames
  
  # Calculate Fst from VC for each pair
  Fst.ls <- lapply(vc, function(pair) {
    
    sum(pair[, 1]) / sum(pair)
    
  })
  
  # Create matrix of pairwise Fst
  Fst <- array(0, c(nsillys, nsillys), dimnames = list(sillyvec, sillyvec))
  
  Fst[lower.tri(Fst)] <- unlist(Fst.ls)
  
  Fst[upper.tri(Fst)] <- t(Fst)[upper.tri(Fst)]
  
  # Create tree
  tree <- ape::nj(Fst)

  # Bootstraps
  message("\nCalculate bootstrap Fst values\n", sep = '')
  
  trees.bootstrapFst <- foreach::foreach(i = seq(nboots), .packages = "ape") %dopar% {
    
    temploci <- sample(loci, nloci, replace = TRUE)
    
    tempFst.ls <- lapply(vc, function(pair) {
      
      sum(pair[temploci, 1]) / sum(pair[temploci, 1:3])
      
    })
    
    tempFst <- array(0, c(nsillys, nsillys), dimnames = list(sillyvec, sillyvec))
    
    tempFst[lower.tri(tempFst)] <- unlist(tempFst.ls)
    
    tempFst[upper.tri(tempFst)] <- t(Fst)[upper.tri(tempFst)]
    
    list(trees = ape::nj(tempFst), bootstrapFst = tempFst)
    
  }

  bootstrapFst <- lapply(trees.bootstrapFst, function(i) {
    
    i[["bootstrapFst"]]
    
  })
  
  trees <- lapply(trees.bootstrapFst, function(i) {
    
    i[["trees"]]
    
  })
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  message("\nBootstrap tree nodes\n", sep = '')
  
  bootstrap <- ape::prop.clades(tree, trees)

  # Create final PairwiseFstTree list object
  if(returnbootstrapFst) {
    
    PairwiseFstTree <- list(tree = tree, bootstrap = bootstrap, PairwiseFst = Fst, vc = vc, BootstrapFst = bootstrapFst)
    
  } else {
    
    PairwiseFstTree <- list(tree = tree, bootstrap = bootstrap, PairwiseFst = Fst, vc = vc)
    
  }
  
  # Save tree before exiting
  dput(x = PairwiseFstTree, file = paste(dir,"\\", length(sillyvec), "Pops", length(loci), "Loci_", "PairwiseFstTree.txt", sep = ""))
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  print(fulltime) 
  
  return(PairwiseFstTree)

}