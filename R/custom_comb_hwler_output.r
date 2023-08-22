#' Custom Combine HWLER Output
#'
#' This function computes summary statistics from HWLER (Hardy-Weinberg and Linkage Equilibrium Sampler) output.
#'
#' @param groupvec A numeric vector indicating the reporting group affiliation for each population in the baseline, where `groupvec` is the same length as the number of populations in the known baseline..
#' @param groupnames A character vector of group names corresponding to the groups in the baseline, where \code{length(groupnames) == max(groupvec)}.
#' @param maindir Directory where the mixture folders are located, where the results for each mixture must be in their own folder with the same name as the mixture.
#' @param mixvec A character vector of mixture silly codes, used to read in the HWLER output.
#' @param ext Extension of the HWLER output file you want to read in and summarize, can either be "RGN" (reporting group estimates) or "BOT" (population estimates). Use "BOT" when resummarizing.
#' @param burn Proportion of iterations to drop from the beginning of each chain. For example, for 40,000 iterations setting \code{burn = 0.5} (default) will drop the first 20,000 iterations.
#' @param alpha Numeric constant specifying credibility intervals; default is 0.1, which gives 90% CIs (i.e., 5% and 95%).
#' @param PosteriorOutput A logical value indicating whether to output the posterior values.
#'
#' @details 
#' This function makes all accommodations for the extra baseline group. 
#' In addition to summary statistics, this function also produces a pdf file of trace plots for each group and mixture.
#' 
#' @return A list of combined HWLER output for each mixture, including summary statistics and posterior output (if requested).
#'        
#' @examples
#' \dontrun{
#' setwd("V:/Analysis/2_Central/Sockeye/Cook Inlet/Missing Baseline Analysis") 
#' 
#' groupvec <- c(1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 
#'   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5, 4, 4, 5, 5, 
#'   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
#'   6, 6, 6, 7, 7, 7, 7)
#' groupnames <- c("A", "B", "C", "D", "E", "F", "G")
#' 
#' mixvec <-  c("Western06early", "Western06late")
#' 
#' maindir <- "HWLER/output"
#' 
#' custom_comb_hwler_output(groupvec = groupvec69, groupnames = groupnames, maindir = maindir, mixvec = mixvec, ext = "BOT", burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
#' }
#'
#' @export
custom_comb_hwler_output <- function(groupvec, groupnames, maindir, mixvec, ext = "BOT", burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE){

  G <- max(groupvec)

  groupvec <- c(groupvec, G+1)

  groupnames <- c(groupnames, "GroupX")

  Results <- vector("list", length(mixvec))
  names(Results) <- mixvec

  pdf(file <- paste(maindir,"\\TracePlotsHWLER", format(Sys.time(), "%b%d%Y_%H_%M_%S"), ".pdf", sep = ""), width = 11, height = 8, family = "Helvetica", pointsize = 20)

  for(mix in mixvec){

    HWLERoutput <- read.table(paste0(maindir, "/", mix, "/", mix, "Chain1", ext, ".", ext))[,-1]

    nits <- nrow(HWLERoutput)

    Burn <- max(1, floor(burn*nits))

    GroupedHWLERoutput <- array(NA, c(nits, G+1), dimnames = list(1:nits, groupnames))

    for(group in groupnames){

      g <- match(group, groupnames)

      if(sum(groupvec==g) > 1){

        GroupedHWLERoutput[1:nits,group] <- apply(HWLERoutput[1:nits, groupvec==g], 1, sum)
     
      }      

      if(sum(groupvec==g) == 1){

        GroupedHWLERoutput[1:nits, group] <- HWLERoutput[1:nits, groupvec==g]
     
      } 
   
      plot(GroupedHWLERoutput[1:nits, group], type = "l", ylim = c(0, 1), yaxp = c(0, 1, 10), col = "red", xlab = "Iteration", ylab = "Proportion", main = paste("Trace of ", group, " from ", mix, sep = ""))

      abline(v = Burn, col = "blue")     

    }

    UsedHWLERoutput <- GroupedHWLERoutput[(Burn + 1):nits, groupnames]

    Results[[mix]] <- vector("list", PosteriorOutput + 1)

    if(PosteriorOutput){
   
      names(Results[[mix]]) <- c("SummaryStats","PosteriorOutput")

      Results[[mix]][["SummaryStats"]] <- data.frame(mean = apply(UsedHWLERoutput, 2, mean), sd = apply(UsedHWLERoutput, 2, sd), LowCI = apply(UsedHWLERoutput, 2, function(clmn){quantile(clmn, alpha/2)}), UpCI = apply(UsedHWLERoutput, 2, function(clmn){quantile(clmn, 1-alpha/2)})) 
   
      Results[[mix]][["PosteriorOutput"]] <- UsedHWLERoutput
      
    } 

    if(!PosteriorOutput){
   
      names(Results[[mix]]) <- c("SummaryStats")

      Results[[mix]][["SummaryStats"]] <- data.frame(mean = apply(UsedHWLERoutput, 2, mean), sd = apply(UsedHWLERoutput, 2, sd), LowCI = apply(UsedHWLERoutput, 2, function(clmn){quantile(clmn,alpha/2)}), UpCI = apply(UsedHWLERoutput, 2, function(clmn){quantile(clmn, 1-alpha/2)})) 
   
    }
   
  }

  dev.off()

  return(Results) 

}