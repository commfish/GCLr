#' Create Rubias Mixture and Baseline Files
#'
#' This function creates rubias mixture and baseline files for different test scenarios. These files are supplied to [GCLr::run_rubias_base_eval()] for producing baseline evaluation estimates.
#' 
#' @param sillyvec a vector of silly codes without the ".gcl" extension (e.g., sillyvec <- c("KQUART06","KQUART08","KQUART10"))
#' 
#' @param group_names A character vector of group names with a length equal to the maximum value in `groupvec`
#' 
#' @param loci A character vector of locus names as they are spelled in LOKI.
#' 
#' @param groupvec A numeric vector indicating the group affiliation of each population in `sillyvec`.
#' 
#' @param sample_sizes A tibble produce by [GCLr::base_eval_sample_sizes()] containing four variables (see details): 
#'    \itemize{
#'        \item `test_group`
#'        \item `scenario`
#'        \item `repunit`
#'        \item `samps`
#'    }
#
#' @param test_groups a character vector of groups to test
#' 
#' @param prprtnl If set to TRUE, the samples for each repunit will be selected in proportion to the number of samples in each population (default: FALSE).
#' 
#' @param base.path The file path where the baseline files will be written.
#' 
#' @param mix.path The file path where the mixture files will be written.
#' 
#' @param seed An integer to set the seed for the random sampler function sample().
#' 
#' @param ncores A numeric vector of length one indicating the number of cores to use.
#' 
#' @param file_type Whether to save the baseline and mixture files as .fst (default and faster) or .csv files.
#'
#' @details
#' This function writes out rubias mixtures and baseline files for each test_group and scenario in sample_sizes. For a given test_group and scenario, the mixture file will contain randomly selected samples for each reporting group as defined in the sample_sizes tibble. The corresponding baseline file will contain all baseline samples except for those selected for the mixture.
#' The function can write .csv or .fst files.  .fst files are compressed, so they save hard drive space, and they are faster to save and read back into R. .csv is also an option to make the function backwards compatible with older analyzes that produced .csv files.
#'
#' @examples
#' \dontrun{
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' group_names <- GCLr::ex_baseline$repunit %>%
#'   unique()
#' 
#' groupvec <- GCLr::ex_baseline %>%
#'   dplyr::group_by(collection) %>%
#'   dplyr::filter(dplyr::row_number()==1) %>%
#'   dplyr::pull(repunit) %>%
#'   factor() %>%
#'   as.numeric()
#' 
#' loci <- GCLr::ex_LocusControl$locusnames
#' 
#' sample_sizes <- GCLr::base_eval_sample_sizes(sillyvec = sillyvec, group_names = group_names, groupvec = groupvec, scenarios = round(seq(.01, 1, .01), 2), mixsize = 200, maxprop = 0.5)
#' 
#' GCLr::create_rubias_base_eval(sillyvec = sillyvec, group_names = group_names, test_groups = group_names, loci = loci, groupvec = groupvec, sample_sizes = sample_sizes, prprtnl = TRUE, seed = 123, ncores = 8)
#' 
#' }
#'
#' @aliases CreateRubiasBaselineEval.GCL
#' 
#' @export
create_rubias_base_eval <- function(sillyvec, group_names, loci, groupvec, sample_sizes, test_groups = group_names, prprtnl = FALSE, base.path = "rubias/baseline", mix.path = "rubias/mixture", seed = 123, ncores = 4, file_type = c("fst", "csv")[1]){
  
  if(sum(is.na(match(loci, LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  start_time <- Sys.time()
  
  if(sum(stringr::str_detect(group_names, "\\W"))>0){
    
    stop("Special characters and spaces were detected in your group_names. 
          Using spaces and delimiters other than underscore in your group names may cause function errors later in your analysis.")
    
  }
  
  group_check <- match(test_groups, sample_sizes$test_group %>% unique()) %>% 
    is.na() %>% 
    sum()
  
  if(group_check > 0){
    
    stop(paste("sample_sizes does not contain scenarios for all test_groups", sep = ""))
    
  }
  
  repunit_check <- match(group_names, sample_sizes$repunit %>% unique()) %>% 
    is.na() %>% 
    sum()
  
  if(repunit_check > 0){
    
    stop(paste("sample_sizes does not contain sample sizes for all group_names supplied", sep = ""))
    
  }
 
  #Write out baseline and mixture files for each test_group scenario in sample_sizes
  
  full_base <- create_rubias_base(sillyvec = sillyvec, loci = loci, group_names = group_names, groupvec = groupvec, file = base.path, baseline_name = "full_base")
    
  # Multicore loop
  cl <- parallel::makePSOCKcluster(ncores)
    
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  doRNG::registerDoRNG(seed = seed, once = TRUE) # This sets the seed for the %dorng% loop.
    
  `%dorng%` <- doRNG::`%dorng%`
  
  foreach::foreach(g = test_groups, .packages = c("tidyverse", "rubias")) %dorng% {
    
    #Test group scenarios
    tgscn <- sample_sizes %>% 
      dplyr::filter(test_group == g) %>% 
      dplyr::pull(scenario) %>% 
      unique()
    
   #Loop through scenarios  
   for(scn in tgscn){
      
      my.sample_sizes <- sample_sizes %>% 
        dplyr::filter(test_group == g, scenario == scn)
      
      mixture <- lapply(my.sample_sizes$repunit, function(repu){
        
        n <- my.sample_sizes %>% 
          dplyr::filter(repunit== repu) %>% 
          dplyr::pull(samps)
        
        IDs <- full_base %>% 
          dplyr::filter(repunit==repu) %>% 
          dplyr::select(indiv, collection)
        
        if(prprtnl){
        
          # This function rounds the number of samples to select so they add up to the original n.
          smart.round <- function(x) {
            
            y <- floor(x)
            
            indices <- tail(order(x-y), round(sum(x)) - sum(y))
            
            y[indices] <- y[indices] + 1
            
            y
            
          }
          
          props <- IDs %>% 
            dplyr::group_by(collection) %>% 
            dplyr::summarize(total_n = length(collection), .groups = "drop_last") %>% 
            dplyr::mutate(prop = total_n/sum(total_n), sample_n = smart.round(!!n*prop))
          
          mix_samples <-  lapply(props$collection, function(col){
            
            IDs %>% 
              dplyr::filter(collection==col) %>% 
              dplyr::pull(indiv) %>% 
              sample(size = props %>% filter(collection == col) %>% pull(sample_n))
            
          }) %>% unlist()
          
        } else{mix_samples <- sample(IDs, size = n)}
        
        full_base %>% 
          dplyr::filter(repunit == repu) %>% 
          dplyr::filter(indiv%in%mix_samples)
        
      }) %>% 
        dplyr::bind_rows() %>% 
        dplyr::mutate(sample_type = "mixture", repunit = NA, collection = paste(g, scn, sep = "_"))
      
      baseline <- full_base %>% 
        dplyr::filter(!indiv%in%mixture$indiv)
      
      if(file_type == "csv"){
        
        readr::write_csv(mixture, file = paste0(mix.path, "/", g, "_", scn, ".mix.csv"))
        
        readr::write_csv(baseline, file = paste0(base.path, "/", g, "_", scn, ".base.csv"))
        
      }
      
      if(file_type == "fst"){
        
        fst::write_fst(mixture, path = paste0(mix.path, "/", g, "_", scn, ".mix.fst"), compress = 100)
        
        fst::write_fst(baseline, path = paste0(base.path, "/", g, "_", scn, ".base.fst"), compress = 100)
        
      }
       
    } #End for loop
    
  } #End multicore
  
  parallel::stopCluster(cl)
  
  Sys.time()-start_time
  
}