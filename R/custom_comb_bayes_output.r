#' Custom Combine BAYES Output
#'
#' This function computes summary statistics from `BAYES` output, similar to custom_comb_rubias_output().
#'
#' @param groupvec Numeric vector indicating the reporting group affiliation.
#'                If summarizing for the original reporting groups using the RGN output, groupvec is seq(length(group_names)).
#'                If summarizing for new groups using the BOT output, groupvec is the length of the number of populations in the baseline.
#' @param group_names Character vector of group names where length(group_names) == max(groupvec).
#' @param maindir The directory where the mixture folders are located.
#'                **Note: the results for each mixture must be in their own folder with the same name as the mixture.**
#' @param mixvec Character vector of mixture sillys, used to read in the BAYES output (.RGN or .BOT).
#' @param ext The extension of the BAYES output file you want to read in and summarize, can either be "RGN" (reporting group estimates) or "BOT" (population estimates).
#'                **Note: use "BOT" when resummarizing.**
#' @param nchains The number of chains to summarize (there should be output files for each chain).
#' @param burn The proportion of iterations to drop from the beginning of each chain.
#'                For example, for 40,000 iterations setting burn = 0.5 will drop the first 20,000 iterations.
#' @param alpha Numeric constant specifying credibility intervals, default is 0.1, which gives 90% CIs (i.e. 5% and 95%).
#' @param threshold Numeric constant specifying how low stock comp is before assuming 0, used for `P=0` calculation.
#' @param plot_trace Logical switch, when TRUE will create a trace plot for each mixture chain and repunit (reporting group). 
#'                The trace plots can take a while to render depending on the number of mixtures, reporting groups, chains, and iterations you are summarizing.
#' @param ncores A numeric vector of length one indicating the number of cores to use.
#'
#' @return Returns a tibble with 8 fields for each mixture and repunit (reporting group):
#'   - mixture_collection: Factor of mixtures (only a factor for ordering, plotting purposes).
#'   - repunit: Factor of reporting groups (only a factor for ordering, plotting purposes).
#'   - mean: Mean stock composition.
#'   - sd: Standard deviation.
#'   - median: Median stock composition.
#'   - loCI: Lower bound of credibility interval.
#'   - hiCI: Upper bound of credibility interval.
#'   - P=0: The proportion of the stock comp distribution that was below `threshold` (i.e., posterior probability that stock comp = 0).
#'   - gr: The Gelman-Ruban shrink factor for checking among-chain convergence.
#' 
#' Also returns a trace plot for each mixture and repunit if `plot_trace = TRUE`.
#'
#' @seealso custom_comb_rubias_output()
#'
#' @examples
#' \dontrun{
#' 
#'   setwd("V:/Analysis/2_Central/Sockeye/Cook Inlet/2012 Baseline/Mixture/2016 UCIfisheryMixtures")
#'   attach("2016UCIfisheryMixtureAnalysis.RData")
#'   groupvec <- groupvec
#'   group_names <- groups
#'   mixvec <- mixvec
#'   detach()
#'   
#'   maindir <- "BAYES/Output"
#'
#'   results <- custom_comb_bayes_output(groupvec = groupvec, group_names = group_names, maindir = maindir, mixvec = mixvec, ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, threshold = 5e-7, plot_trace = TRUE, ncores = 8)
#' 
#' }
#'
#' @export
custom_comb_bayes_output <- function(groupvec, group_names, maindir, mixvec, ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, threshhold = 5e-7, plot_trace = FALSE, ncores = 4){
  
  #Error checking
  
  if(length(group_names)!= max(groupvec)){
    
    stop(paste0(length(group_names), " group_names were supplied but there are ", max(groupvec)," groups in the groupvec supplied. Make sure that length(group_names)==max(groupvec)."))
  
    }
  
  G <- max(groupvec)

  C <- length(groupvec)

  nummix <- length(mixvec)
  
  mix_files <- lapply(mixvec, function(mix){ #mixvec loop
    
    paste0(maindir, "/", mix, "/", mix, "Chain", 1:nchains, ext, ".", ext)
    
  }) %>% set_names(mixvec)
  
  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  `%dopar%` <- foreach::`%dopar%`
    
  output0 <- foreach::foreach(mix = mixvec, .packages = c("tidyverse")) %dopar% { #mixvec loop
      
    filenames <- mix_files[[mix]]
 
    files0 <- lapply(filenames, function(filename){
      
      file0 <- as.matrix(read.table(filename)[,-1]) 
      
      if(dim(file0)[[2]] != length(groupvec)){
        
        stop("The length of the groupvec must equal the number of columns in the ", ext, " file: ", dim(file0)[[2]])
        
      }
      
      file0 %>% 
        t() %>% 
        rowsum(group = groupvec) %>% 
        t()
        
      }) 
    
    iter_check <- suppressMessages((lapply(files0, FUN = dim) %>% 
                                      dplyr::bind_cols())[1,] %>% t())
   
    if(length(unique(iter_check[, 1]))>1){
      
      stop("Chains must be the same length!!!")
      
    }
    
    nitter <- dim(files0[[1]])[1]
    
    files <- lapply(files0, function(x){
      
      x %>% 
        tibble::as_tibble() %>% 
        dplyr::mutate(iteration = seq(nitter)) %>% 
        dplyr::filter(iteration > nitter*burn)
      
    })
    
    files4GR <- lapply(files, function(file) {

      coda::as.mcmc(file %>% select(-iteration))

    })%>% coda::as.mcmc.list()

    if(nchains<2){

      GR <- rep(NA, G)

    }else{    

      GR <- coda::gelman.diag(files4GR, multivariate = FALSE, transform = TRUE)[[1]][,1]

    }
    
    posterior <- Reduce(rbind, lapply(files4GR, as.matrix))
    
    trace_out <- lapply(files0, function(file){
      
      file %>% 
        tibble::as_tibble() %>% 
        dplyr::mutate(iteration = 1:nitter)
      
      }) %>% 
      purrr::set_names(1:nchains) %>% 
      bind_rows(.id = "Chain") %>% 
      magrittr::set_colnames(c("Chain", group_names, "iteration"))
    
    list(output = tibble(mixture_collection = mix, 
           repunit = group_names,
           mean = apply(posterior, 2, mean), 
           sd = apply(posterior, 2, sd), 
           median = apply(posterior, 2, median), 
           loCI = apply(posterior, 2, quantile, probs = alpha/2), 
           hiCI = apply(posterior, 2, quantile, probs = 1-alpha/2), 
           `P=0`= apply(posterior, 2, function(clm){sum(clm<threshhold)/length(clm)}), 
           gr = GR),
         trace = trace_out)
    
    } %>% purrr::set_names(mixvec)#End mixvec loop
  
  parallel::stopCluster(cl)

  output <- lapply(1:length(mixvec), function(mix){
    
    output0[[mix]][["output"]]
    
    }) %>% bind_rows()
  
  if(plot_trace){
    
    trace <- lapply(1:length(mixvec), function(mix){
      
      output0[[mix]][["trace"]] 
      
    }) %>% 
      purrr::set_names(mixvec) %>%
      dplyr::bind_rows(.id = "mixture_collection") %>% 
      tidyr::pivot_longer(all_of(group_names), names_to = "repunit", values_to = "est") %>% 
      dplyr::mutate(Chain = paste0("Chain ", Chain))
    
    # Plot each mixture separately and facet by repunit and chain 
    # These can take quite a while to plot depending on the number of mixtures, iterations, and groups.
    
    for(mix in mixvec){ 
      
      my.trace <- trace %>% 
        dplyr::filter(mixture_collection == mix) %>% 
        dplyr::mutate(repunit = factor(x = repunit, levels = group_names))
      
       canvas <- my.trace %>% 
              dplyr::filter(mixture_collection == mix) %>% 
              dplyr::mutate(repunit = factor(x = repunit, levels = group_names)) %>%  # order repunit
              ggplot2::ggplot(aes(x = iteration, y = est, colour = repunit)) 
        
      print(canvas +
         ggplot2::geom_line() +
         ggplot2::ylim(0, 1) +
         ggplot2::geom_vline(xintercept = burn*max(my.trace$iteration)) +
         ggplot2::theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
         ggplot2::facet_grid(repunit ~ Chain)+
         ggplot2::ggtitle(label = mix)
      )
              
    }
    
  }

  return(output)

}

