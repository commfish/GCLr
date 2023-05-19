create_rubias_base_eval <- function(sillyvec, group_names, loci, groupvec, sample_sizes, test_groups = group_names, prprtnl = FALSE, base.path = "rubias/baseline", mix.path = "rubias/mixture", seed = 123, ncores = 4, file_type = c("fst", "csv")[1]){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates rubias mixture and baseline files for different proof test scenarios.  
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   group_names - character vector of group names the length of max(groupvec)
  #
  #   groupvec - numeric vector indicating the group affiliation of each pop in sillyvec
  #
  #   loci - character vector of locus names as they are spelled in LOKI - example: c("GTH2B-550", "NOD1", "Ots_100884-287")
  #
  #   test_groups - a character vector of groups to test
  #
  #   sample_sizes - a tibble containing 4 variables test_group, scenario, repunit, and samps
  #    
  #    Here's an example tibble for one mixture with 1% of samples going to the Upper_Susitna repunit (test_group) and the remaining  
  #    samples spread among the 5 other repunits.
  #     test_group    scenario repunit         samps
  #     <chr>            <dbl> <chr>           <dbl>
  #     Upper_Susitna     0.01 Upper_Susitna       2  #Notice Upper_Susitna has 2 samples selected for a 200 sample mixuture (2/200 = 0.01)
  #     Upper_Susitna     0.01 Chulitna           39  #The remaining 198 samples are spread across the other reporting groups (repunits) 
  #     Upper_Susitna     0.01 Talkeetna          47
  #     Upper_Susitna     0.01 Eastern_Susitna    40
  #     Upper_Susitna     0.01 Deshka             35
  #     Upper_Susitna     0.01 Yentna             37
  #
  #    prprtnl - logical, if set to TRUE the samples for each repunit will be selected in proportion to the number of samples in each population. 
  #              Setting prprtnl = TRUE helps avoid oversampling populations when a reporting group contains only a few populations with small sample sizes.
  #
  #    base.path - the file path where the baseline .csv files will be written
  #
  #    mix.path - the file path where the mixture .csv files will be written
  #
  #    seed - integer to set the seed for the random sampler function sample()
  #
  #    ncores - a numeric vector of length one indicating the number of cores to use
  #   
  #    file_type - whether you want the baseline and mixture files saved as .fst (default and much faster!) or .csv files. 
  #                This was argument was added for backwards compatibility.  
  #
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function writes out rubias mixtures and baseline .csv files for each test_group and scenario in sample_sizes. 
  #   For a given test_group and scenario, the mixture file will contain randomly selected samples for each reporting group as defined in the sample_sizes tibble. 
  #   The corresponding baseline file will contain all baseline samples except for those selected for the mixture.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  attach("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  #  Final_Pops <- Final_Pops %>% mutate(group = factor(group, levels = unique(group)))
  #  sample_sizes <- base_eval_sample_sizes(sillyvec = Final_Pops$silly, group_names = Final_Pops$group %>% levels(), groupvec = Final_Pops$group %>% as.numeric(), scenarios = round(seq(.01, 1, .01), 2), mixsize = 200, maxprop = 0.5)
  # 
  #  create_rubias_base_eval(sillyvec = Final_Pops$silly, group_names = Final_Pops$group %>% levels(), test_groups = (Final_Pops$group %>% levels())[1:2], loci = loci80, groupvec = Final_Pops$group %>% as.numeric(), sample_sizes = sample_sizes, prprtnl = TRUE, seed = 123, ncores = 8)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  if(sum(is.na(match(loci, LocusControl$locusnames)))){
    
    stop(paste("'", loci[is.na(match(loci, LocusControl$locusnames))], "' from argument 'loci' not found in 'LocusControl' object!!!", sep = ""))
    
  }
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  start_time <- Sys.time()
  
  if(sum(str_detect(group_names, "\\W"))>0){
    
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