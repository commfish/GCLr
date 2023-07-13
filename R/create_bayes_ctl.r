#' @title Create Bayes Control Files
#'
#' @description
#' This function creates BAYES control (.ctl) files for each mixture in `mixvec` and each chain in `nchains`. It takes
#' various arguments such as the population sillys used to create the baseline file, the loci used to produce the baseline
#' and mixture files, the mixtures for which control files are to be created, the directory where the control files will
#' be saved, the baseline name (without the .bse extension), the number of MCMC iterations, the number of MCMC chains, the
#' group numbers and priors for each silly, the initial starting values for each chain, the random seeds, the thinning
#' intervals for MCMC sampling, the Fortran format of the mixture and baseline files, and various logical switches.
#'
#' @param sillyvec Character vector of population sillys used to create the baseline file.
#' @param loci Character vector of the loci used to produce the baseline and mixture files.
#' @param mixvec Character vector of mixture ".gcl" objects for which control files are to be created.
#' @param dir Character vector of the directory where the control files will be saved.
#' @param baseline_name Character string of the baseline file name without the .bse extension.
#' @param nreps Numeric; the number of MCMC iterations for the mixtures.
#' @param nchains The number of MCMC chains for the mixtures.
#' @param groupvec Numeric vector of group numbers corresponding to each silly in `sillyvec`.
#' @param priorvec Numeric vector of priors corresponding to each silly in `sillyvec`.
#' @param initmat Numeric matrix of initial starting values for each chain with, nrow(`initmat`) == length(`sillyvec`) & ncol(`initmat`) == nchains.
#' @param seeds Matrix of random seeds containing 3 seeds per chain with, nrow(`seeds`) == 3 & ncol(`seeds`) == `nchains`.
#' @param thin Thinning intervals for MCMC sampling of 1) stock proportions, 2) baseline allele or type relative frequencies, and 3) stock assignments of each mixture individual.
#' @param mixfortran The Fortran format of the mixture files, as created by [GCL::create_bayes_mix()].
#' @param basefortran The Fortran format of the baseline file,as created by [GCL::create_bayes_base()].
#' @param switches Character string of logical switches:
#'   - Switch 1: Baseline file printed (default: "F")
#'   - Switch 2: Mixture file printed (default: "T")
#'   - Switch 3: MCMC sample of baseline relative frequencies printed (default: "F")
#'   - Switch 4: Baseline file contains counts, not relative frequencies (default: "T")
#'   - Switch 5: Stock assignment distribution for each individual in the mixture sample printed (default: "F")
#'   - Switch 6: Stock-group or regional estimates printed (default: "T")
#'   - Switch 7: Individual stock estimates are suppressed (default: "F")
#'   
#' @returns This function writes out bayes control (.ctl) files to `dir`
#' 
#' @examples
#' \dontrun{
#' sillyvec <- c("KSUSC18FW", "KSUSCN18", "KSUSC19FW", "KSUSCN19")
#' loci <- c("locus1", "locus2", "locus3")
#' mixvec <- c("mixture1", "mixture2")
#' baseline_name <- "baseline"
#' nreps <- 40000
#' nchains <- 5
#' groupvec <- c(1, 2, 3, 4)
#' priorvec <- c(0.1, 0.2, 0.3, 0.4)
#' initmat <- matrix(0.1, nrow = 4, ncol = 5)
#' dir <- "control_files"
#' seeds <- matrix(12345, nrow = 3, ncol = 5)
#' thin <- c(1, 1, 1)
#' mixfortran <- "(1X,I6)"
#' basefortran <- "(1X,I6)"
#' switches <- "F T F T F T F"
#' create_bayes_ctl(sillyvec, loci, mixvec, baseline_name, nreps, nchains, groupvec, priorvec, initmat, dir, seeds, thin, mixfortran, basefortran, switches)
#' }
#' 
#' @seealso See bayes manual for addtional details:  [V:\Software\BAYES\MANUAL.DOC]
#' @export
create_bayes_ctl <- function(sillyvec, loci, mixvec, baseline_name, nreps = 40000, nchains, groupvec, priorvec, initmat, dir, seeds = matrix(sample(seq(10000), 3 * nchains), nrow = 3), thin = c(1, 1, 1), mixfortran, basefortran, switches = "F T F T F T F") {

  if (!all(loci %in% LocusControl$locusnames)) {
    
    stop(paste0("\n'", setdiff(loci, LocusControl$locusnames), "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
  }
  
  sapply(mixvec, function(mix){
    
    chains <- paste("Chain", 1:nchains, sep = "")
    
    dirs <- paste(dir, "/", mix, "Chain", 1:nchains, ".ctl", sep = "") %>% setNames(chains)

    priorvec <- cbind(substring(format(round(priorvec, 6), nsmall = 6), first = 2))
    
    initmat <- cbind(substring(format(round(initmat, 6), nsmall = 6), first = 2))
    
    if (nchains > 1) {
      
      dimnames(initmat)[[2]] <- chains
      
    }
    
    nsillys <- length(sillyvec)
    
    ploidy <- LocusControl$ploidy[loci]
    
    nalleles <- LocusControl$nalleles[loci]
    
    nloci <- length(loci)
    
    seeds <- cbind(seeds)
    dimnames(seeds) <- list(1:3, chains)
    
    files <- lapply(chains, function(chain){paste0(mix, chain)}) %>% setNames(chains)
      
    files <- lapply(files, function(file){rbind(file, paste0(baseline_name, ".bse"))}) %>% setNames(chains)
  
    files <- lapply(files, function(file){rbind(file, paste0(mix, ".mix"))}) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){rbind(files[[chain]], paste0(mix, chain, "SUM.SUM"))}) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){rbind(files[[chain]], paste0(mix, chain, "BOT.BOT"))}) %>% setNames(chains)
 
    files <- lapply(chains, function(chain){rbind(files[[chain]], paste0(mix, chain, "FRQ.FRQ"))}) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){rbind(files[[chain]], paste0(mix, chain,"BO1.BO1"))}) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){rbind(files[[chain]], paste0(mix, chain, "CLS.CLS"))}) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){rbind(files[[chain]], paste0(mix, chain, "RGN.RGN"))}) %>% setNames(chains)
     
    files <- lapply(chains, function(chain){rbind(files[[chain]], nreps)}) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){rbind(files[[chain]], length(sillyvec))}) %>% setNames(chains)
 
    files <- lapply(chains, function(chain){rbind(files[[chain]], length(loci))}) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){rbind(files[[chain]], cbind(seeds[, chain]))}) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){rbind(files[[chain]], cbind(thin))}) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){rbind(files[[chain]], mixfortran)}) %>% setNames(chains)
 
    files <- lapply(chains, function(chain){rbind(files[[chain]], basefortran)}) %>% setNames(chains)

    files <- lapply(chains, function(chain){rbind(files[[chain]], switches)}) %>% setNames(chains)
    
    files <- lapply(seq(length(chains)), function(chain){
      
      rbind(files[[chain]], cbind(sapply(1:nloci, function(d){paste(sprintf("%3s", d), sprintf("%2s", nalleles[loci[d]]), sprintf("%2s", ifelse(ploidy[loci[d]] == 2, "T", "F")), "", loci[d], collapse = "")})))
      
      }) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){
      
      rbind(files[[chain]], cbind(sapply(1:nsillys, function(i){
        
        paste(sprintf("%3s", i), sprintf("%2s", groupvec[i]), sprintf("%7s", priorvec[i]), "", format(ifelse(nchar(sillyvec[i]) < 18, sillyvec[i], substr(sillyvec[i], start = 1, stop = 18)), width = 18), sprintf("%7s", initmat[i, chain]), collapse = "")
        
        })))
      
      }) %>% setNames(chains)

    empty <- sapply(chains, function(chain){write.table(files[[chain]], dirs[chain], quote = FALSE, row.names = FALSE, col.names = FALSE)})
 
    }) # End mix loop   
  
}