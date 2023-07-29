#' Get Difference in Allele Frequencies Between Baseline and BAYES Posteriors
#' 
#' This function reads in BAYES MCMC samples of baseline frequencies for each population (one file per chain) and compares them to the original baseline allele frequencies (see details).
#' 
#' @param sillyvec sillyvec Character vector of population sillys used to create the baseline file.
#' 
#' @param loci Character vector of the loci used to produce the baseline and mixture files.
#' 
#' @param mydir the director where the BAYES mixture output files are located.  e.g., "BAYES/output" This directory should contain a folder named `mixname`, where the .FRQ files are saved.
#' 
#' @param nchains the number of MCMC chains analyzed for the mixture
#' 
#' @param mixname the name of the mixture. There should be a folder in `mydir` named mix
#' 
#' @details 
#' During a mixed stock analysis in BAYES, the original allele frequencies get updated at each MCMC iteration. 
#' Looking at the difference between the average MCMC allele frequency and the original baseline frequencies can be useful when looking into reasons for non-convergence among chains (i.e., Gelman-Ruban shrink factor > 1.2) 
#' Populations/reporting groups with large changes in allele frequencies for certain population/reporting groups can indicate that the mixture contains fish from a distinct population that is not represented in the baseline.
#' 
#' @seealso See [BAYES manual](system.file("BAYES", "MANUAL.DOC", package = "GCLr")) for additional details.
#' 
#' @returns a list containing 4 elements:
#'   \itemize{
#'     \item \code{Outliers}: a logical matrix with `nrows = length(sillyvec)` and `ncols = nchains` indicating where a population was an outlier for a given chain.
#'     \item \code{MeanAbsDiffByPop}: a list with one element per chain containing a vector of mean (among loci) absolute allele frequency differences for each population in `sillyvec`
#'     \item \code{PriorAlleleFreqs}: a 3 dimensional array of baseline allele frequencies for each locus in `loci`, silly in `sillyvec`, and allele.
#'     \item \code{PosteriorAlleleFreqs}: a list with one element per chain containing a vector of BAYES MCMC allele frequencies for each population in `sillyvec` and each locus in `loci`
#'   }
#'   
#' list(Outliers = Outliers, MeanAbsDiffByPop = myabsdiffsByPop, PriorAlleleFreqs = rawQ, PosteriorAlleleFreqs = myQ)
#' @examples
#' \dontrun{
#' 
#' post_prior_diff_allele_freq(sillyvec = sillyvec, loci = loci, mydir = "BAYES/output", nchains = 5, mixname = "my.mixname")
#' 
#' }
#' 
#' @export
post_prior_diff_allele_freq <- function(sillyvec, loci, mydir, nchains, mixname){

  require("outliers")

  baseline <- calc_freq_pop(sillyvec,loci)

  loci <- dimnames(baseline[1, , ])[[1]]

  nalleles <- apply(baseline[1, , ], 1, function(vec){sum(!is.na(vec))})

  C <- dim(baseline)[1]

  L <- dim(baseline)[2]

  rawQ <- baseline

  for(i in 1:C){
    
    for(loc in 1:L){
      
       rawQ[i, loc, 1:nalleles[loc]] <- (baseline[i, loc, 1:nalleles[loc]] + 1/nalleles[loc])/sum(baseline[i, loc, 1:nalleles[loc]] + 1/nalleles[loc]) 
       
    }
    
  }

  ChainNames <- paste0("Chain", 1:nchains)

  mydirs <- paste(mydir, "\\", mixname, ChainNames, "FRQ.FRQ", sep = "")

  myfiles <- lapply(mydirs, function(dr){
    
    myfile <- readLines(dr); 
    
    myfile[(length(myfile)/2+1):length(myfile)]
    
    })

  names(myfiles) <- ChainNames

  myQ <- lapply(myfiles,function(myfile){
    
        ans <- t(vapply(myfile, function(strng){
              vec <- strsplit(strng," ")[[1]];
              vec <- vec[vec!=""];
              vec2 <- rep(NA,max(nalleles)+4);
              vec2[1:length(vec)] <- vec;
              as.numeric(vec2)
            }, as.numeric(rep(NA, 4+max(nalleles)))));
        
        dimnames(ans) <- list(1:nrow(ans), c("Iteration", "Pop", "Locus", "n", paste("Allele", 1:max(nalleles), sep = "")));
        
        ans <- data.frame(ans);
        
        q <- rawQ;
        
        for(i in 1:C){
          
          for(loc in 1:L){
            
            q[i, loc, 1:nalleles[loc]] <- apply(ans[ans$Pop==i & ans$Locus==loc, 5:(4+nalleles[loc])], 2, mean)
            
          }
          
        };
        
        q;
        
        })

  names(myQ) <- ChainNames

  myabsdiffsByPopByLocus <- sapply(myQ,function(qq){ 
    
    sapply(loci, function(locus){apply(abs(qq[1:C, locus, 1:nalleles[locus]]-rawQ[1:C, locus, 1:nalleles[locus]]), 1, mean) 
      
      }) 
    
    },simplify=FALSE)

  myabsdiffsByPop <- sapply(myabsdiffsByPopByLocus, function(diffs){
    
    apply(diffs,1,mean)
    
    },simplify=FALSE)

  Outliers <- sapply(myabsdiffsByPop, function(popdiffs){
    
    outliers::outlier(popdiffs, logical = TRUE)
    
    })

  return(list(Outliers = Outliers, MeanAbsDiffByPop = myabsdiffsByPop, PriorAlleleFreqs = rawQ, PosteriorAlleleFreqs = myQ))
  
}