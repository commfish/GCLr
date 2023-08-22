#' Create HWLER Control Files
#'
#' This function creates HWLER (Hardy-Weinberg and Linkage Equilibrium Sampler) control files for each chain based on the specified parameters.
#'
#' @param sillyvec a character vector of baseline sillys included in the input file. Leave this `NULL` if no baseline individuals are present in the input file.
#' @param loci  character vector of the loci used to produce the baseline and mixture files.
#' @param input The mixture and/or baseline input file name without the extension. This file can contain mixture data only, baseline data only, or a combination mixture and baseline. (see [HWLER manual](system.file("HWLER", "HWLER_manual.doc", package = "GCLr")) for more details) 
#' @param nsamples a vector specifying the number of samples for each chain
#' @param nchains the number of MCMC chains to analyze the mixtures
#' @param dir the directory path where the control files will be saved
#' @param initval the initial starting value for each MCMC chain 
#' @param seeds a matrix of random seeds containing 3 seeds per chain, where \code{nrow(seeds) = 3} and \code{ncol(seeds) = nchains}.
#' @param thin thinning intervals for MCMC sample of 1) stock proportions, 2) baseline allele or type relative frequencies, and 3) stock assignments of each mixture individual.
#' @param inputfortran the FORTRAN format of the input created by [GCLr::create_hwler_input()]
#' @param switches a character string of logical switches, with default value "T T T F T T T T T" (see details).
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()], (default = LocusControl)  
#' 
#' @details 
#' The HWLER input file can contain only mixture individuals, only baseline individuals, or a combination of mixture and baseline individuals depending on the type of analysis you are doing. 
#' The `switches` argument has 9 program options turned “on” with “T” for true or “off” with “F” for false:
#'   \enumerate{
#'     \item Baseline printed
#'     \item Mixture printed
#'     \item Output stock assignments of individuals (baseline & mixture) or snapshots of cluster assignments of individuals (baseline or mixture only)
#'     \item Suppress the printing of MCMC samples of stock composition (baseline & mixture) or snapshots of numbers of individuals assigned to clusters (baseline or mixture only)
#'     \item Baseline provided 
#'     \item Mixture provided 
#'     \item Individual identification provided in data file
#'     \item Output cluster information for the binary tree
#'     \item Sample the Dirichlet mass parameter
#'     }
#'     
#' @seealso See HWLER manual for additional details: [HWLER manual](system.file("HWLER", "HWLER_manual.doc", package = "GCLr"))
#'     
#' @returns Writes out HWLER control (.ctl) files.
#'
#' @examples
#' \dontrun{
#' load("V:/Analysis/2_Central/Sockeye/Cook Inlet/Missing Baseline Analysis/Missing Baseline Analysis.RData")
#' HWLERControlFile.GCL(sillyvec = PopNames69, loci = loci34, input = "myHWLER.input", nsamples = c(22000, 2, 5), nchains = 1, dir = "HWLER/control", initval = "1.0", seeds = matrix(sample(seq(2000000000), 3*nchains), nrow = 3), thin = c(1, 100, 4), inputfortran = inputfortran, switches = "F F T F T T T T F")
#'}
#' @export
create_hwler_ctl <- function(sillyvec = NULL, loci, input, nsamples = c(1000, 2, 5), nchains = 5, dir, initval = "1.0", seeds = matrix(sample(seq(10000), 3*nchains), nrow = 3), thin = c(1, 1, 1), inputfortran, switches = "T T T F T T T T T", LocusCtl = LocusControl){

  if(sum(is.na(match(loci, LocCtl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocCtl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  chains <- paste("Chain", 1:nchains, sep = "")

  dirs <- paste(dir, "\\", input, "Chain", 1:nchains, ".ctl", sep = "")
  names(dirs) <- chains

  nsillys <- length(sillyvec)

  nalleles <- LocusCtl$nalleles[loci]

  nloci <- length(loci)
  
  seeds <- cbind(seeds)
  dimnames(seeds) <- list(1:3, chains)

  file <- lapply(chains, function(chain){paste(input, chain, sep = "")})
  names(file) <- chains

  file <- lapply(file,function(file){rbind(file, paste(input, "input", sep = "."))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], paste(input, chain, "SUM.SUM", sep = ""))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], paste(input, chain, "BOT.BOT", sep = ""))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], paste(input, chain, "CLS.CLS", sep = ""))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], paste(input, chain, "STK.STK", sep = ""))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], paste(input, chain, "TXT.TXT", sep = ""))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], paste(input, chain, "ALP.ALP", sep = ""))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], paste(nsamples[1], nsamples[2], nsamples[3], colapse = "    "))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], length(sillyvec))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], length(loci))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], initval)})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], cbind(seeds[, chain]))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], cbind(thin))})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], inputfortran)})
  names(file) <- chains

  file <- lapply(chains, function(chain){rbind(file[[chain]], switches)})
  names(file) <- chains

  file <- lapply(seq(length(chains)), function(chain){rbind(file[[chain]], cbind(sapply(1:nloci, function(d){paste(sprintf("%3s", d), sprintf("%2s", nalleles[loci[d]]), "", loci[d], collapse = "")})))})
  names(file) <- chains

  if(!is.null(sillyvec)){
    
    file <- lapply(chains, function(chain){rbind(file[[chain]], cbind(sapply(1:nsillys, function(i){paste(sprintf("%3s", i), "", format(ifelse(nchar(sillyvec[i]) < 18, sillyvec[i], substr(sillyvec[i], start = 1, stop = 18)), width = 18), collapse = "")})))})
    names(file) <- chains
    
  }

  empty <- sapply(chains, function(chain){write.table(file[[chain]], dirs[chain], quote = FALSE, row.names = FALSE, col.names = FALSE)})

}