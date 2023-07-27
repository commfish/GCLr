#' Compare Compositions Between Mixtures
#'
#' This function compares the difference in stock compositions between two mixtures.
#' **This needs a better explanation of what it does**
#'
#' @param mixnames A character vector of length 2 containing the mixture names to compare.
#' @param groupnames A character vector giving the names of the stock groups resolved by the mixtures, with length equal to the number columns in the RGN files.
#' @param mixdir An atomic string giving the folder path where the mixture sub-folders are located.
#' @param diffs Differences to check.
#' @param nchains The number of MCMC chains for summarizing the mixture results; default is 5.
#' @param burn The proportion of iterations to burn prior to summarizing results; default is 1/2.
#'
#' @return A list with two elements: \code{one.sided} and \code{two.sided}.
#' \describe{
#'   \item{\code{one.sided}}{A tibble with one-sided p-values for each group.}
#'   \item{\code{two.sided}}{A tibble with two-sided p-values for each group at different \code{diffs}.}
#' }
#'
#' @examples
#'
#' \dontrun{
#' load("V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Mixture/2013 UCIfisheryMixtures/2013UCIfisheryMixtureAnalysis.RData")
#' mixdir <- "V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Mixture/2013 UCIfisheryMixtures/BAYES/Output"
#' groups <- c("group1", "group2", "group3")  # Replace with actual group names
#' compare_comps_between(mixnames = c("DriftExpCorr.Jul11", "Drift.Jul8"), groupnames = groups, mixdir = mixdir, diffs = seq(10)/100, nchains = 5, burn = 1/2)
#' }
#'
#' @export
#' 
compare_comps_between <- function(mixnames, groupnames, mixdir, diffs = seq(10)/100, nchains = 5, burn = 1/2){
  
  if(length(mixnames)>2|length(mixnames)<2){
    
    stop("This function compares estimates between two mixtures. Make sure mixnames contains exactly two mixture names and try again.")
    
  }
  
  dirs <- paste0(mixdir, "/", mixnames)

  chains <- paste0("Chain", seq(nchains), "RGN.RGN")
  
  fnames <- lapply(dirs, function(dir){
    
    files <- list.files(path = dir, pattern = "*.RGN", full.names = TRUE)
      
    if(length(files)!=nchains){
      
     stop("The number of .RGN files does not equal nchains") 
      
    }
    
    files
    
    }) %>% setNames(mixnames)

  files <- sapply(mixnames, function(mix){
    
    lapply(fnames[[mix]], function(file){read.table(file)})
    
    }, simplify = FALSE)
  
  # Remove burnin and collapse down chain outputs into a single posterior distribution for each mixture

  output <- sapply(mixnames, function(mix){ 
    
    Reduce(rbind, lapply(files[[mix]], function(file){
      
      file[seq(floor(nrow(file) * burn) + 1, nrow(file)), -1]
      
      }))
    
    }, simplify = FALSE)

  dsim0 <- Reduce("-", output) #Gets differences in posterior estimates

  colnames(dsim0) <- groupnames
  
  # Create plots of posterior differences by group
  dsim_df <- as_tibble(dsim0) %>% 
    pivot_longer(everything(), names_to = "group", values_to = "difference") %>% 
    mutate(group = factor(group, levels = groupnames))
  
  suppressMessages( print(dsim_df %>% 
    ggplot2::ggplot(aes(x = difference)) + 
    geom_freqpoly() + 
    facet_wrap(~group)+
    geom_vline(xintercept = 0, lty = 2)))

  dsim <- t(diag((-1) ^ (apply(dsim0, 2, mean) < 0)) %*% t(dsim0))

  one.sided <- cbind(one.sided = c(setNames(apply(dsim < 0, 2, mean), groupnames), Overall = mean(apply(dsim < 0, 1, all))), diff.mean = c(apply(dsim0, 2, mean), NA)) %>% 
    as_tibble(rownames = "group")

  two.sided <- sapply(setNames(diffs, diffs), function(dd){c(apply(abs(dsim0) < dd, 2, mean), Overall = mean(apply(abs(dsim0) < dd, 1, all)))}) %>% 
    as_tibble(rownames = "group")  

  return(list(one.sided = one.sided, two.sided = two.sided))

}  



