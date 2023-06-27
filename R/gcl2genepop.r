gcl2genepop <- function(sillyvec, loci, path, VialNums = TRUE, usat = FALSE, ncores = 4, npops = NULL){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Write out a GENEPOP input file from "*.gcl" objects.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   sillyvec - character vector of SILLYs - example: c(KCRES10", "KSTRA10", "KNIKOL12", "KNIKOL13")
  #
  #   loci - character vector of locus names as they are spelled in LOKI - example: c("GTH2B-550", "NOD1", "Ots_100884-287")
  #
  #   path - full file path to write out the GENEPOP file with "\\" or "/" separator between folders 
  #           example: "V:\\Analysis\\2_Central\\Chinook\\Cook Inlet\\2019\\2019_UCI_Chinook_baseline_hap_data\\GENEPOP\\genepopfile.gen"
  #                   or "V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/GENEPOP/genepopfile.gen"
  #           
  #          When npops is supplied, only supply the path to the folder where the files will be written, files will be given generic names.
  # 
  #   VialNums - logical; if TRUE (default), vial numbers will be included for each individual next to their silly code separated by an underscore (e.g. KCRESC10_1)
  #                       if FALSE, only the silly code will be included for each individual.
  #
  #   usat - logical; whether the data are from microsatellites (TRUE) or not (FALSE). 
  #          This is included because GENEPOP only accepts numberic alleles; SNP alleles have to be converted from character to numeric and microsatellite alleles are already numeric.
  #
  #   ncores - the number of cores for mulitcoring using doParallel and foreach. 
  #
  #   npops - optional; the number of pops for each file.  If supplied, multiple files with generic file names will be written out with npops per file. 
  #           For example, if there are 10 pops in sillyvec and npops = 2, the following files will be written to the folder supplied in the path argument: Pops1to2.gen.txt, Pops3to4.gen.txt, Pops5to6.gen.txt, Pops7to8.gen.txt, Pops9to10.gen.txt 
  #           This argument is useful when running Linkage Disequilibrium tests in genepop using the gcltest_LD function so things finish faster. 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Writes out a GENEPOP file to path.
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #  
  #  source(paste0(path.expand("~/R/"), "Functions.R"))#GCL functions
  #
  #  ### SNP
  #  attach("V:/Analysis/2_Central/Coho/Cook Inlet/2019/2019_Cook_Inlet_coho_baseline/2019_Cook_Inlet_coho_baseline.RData")
  #  old2new_locuscontrol()
  #  sapply(sillyvec104, function(silly){assign(paste0(silly, ".gcl"), get(paste0(silly, ".gcl")), pos = -1, envir = .GlobalEnv)})
  #  sillyvec <- sillyvec104
  #  loci <- loci81
  #  detach()
  #  old2new_gcl(sillyvec)
  #  gcl2genepop(sillyvec = sillyvec, loci = loci, path = "GENEPOP/genepopfile.gen", VialNums = TRUE, usat = FALSE, ncores = 8, npops = NULL)
  #
  #  ### uSAT
  #  attach("V:/Analysis/2_Central/Chinook/Susitna-Watana/2017/Upper and Middle River Analysis/SuWuUpperMiddleRiverAnalysis_tidy.RData")
  #  sapply(Sillys, function(silly){assign(paste0(silly, ".gcl"), get(paste0(silly, ".gcl")), pos = -1, envir = .GlobalEnv)})
  #  sillyvec <- Sillys
  #  loci <- loci13
  #  LocusControl <- LocusControl
  #  detach()
  #  gcl2genepop(sillyvec = sillyvec, loci = loci, path = "V:/Analysis/2_Central/Chinook/Susitna-Watana/2017/Upper and Middle River Analysis/Genepop/genepopfile.gen", VialNums = TRUE, usat = TRUE, ncores = 4)
  #
  #  ### SNP with npops argument
  #  attach("C:/2019 coho baseline/2019_Cook_Inlet_coho_baseline_new.RData")
  #  LocusControl <- LocusControl
  #  sapply(sillyvec104, function(silly){assign(paste0(silly, ".gcl"), get(paste0(silly, ".gcl")), pos = -1, envir = .GlobalEnv)})
  #  sillyvec <- sillyvec104
  #  loci <- loci81
  #  detach()
  #  gcl2genepop(sillyvec = sillyvec[1:36], loci = loci, path = "GENEPOP", VialNums = TRUE, usat = FALSE, ncores = 8, npops = 2)
  # 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  start_time <- Sys.time()
  
  if(!exists("LocusControl")){
    
    stop("'LocusControl' not yet built.")
    
  }
  
  if(!all(loci %in% LocusControl$locusnames)){
    
    stop(paste0("'", setdiff(loci, LocusControl$locusnames), "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
  }
  
  if(!is.null(npops) & !dir.exists(path)){
    
    stop("The path supplied is not a folder in the current working directory.\nPlease supply a path to an existing folder when using the npops argument.")
    
  }
  
  if(is.null(npops) & dir.exists(path)){
    
    stop("The path supplied is an existing folder, not a file path. Make sure to supply the full file path. e.g., 'GENEPOP/genepopexample.gen.txt'")
    
  }
  
  
  if(ncores > parallel::detectCores()) {
    
    stop("'ncores' is greater than the number of cores available on machine\nUse 'detectCores()' to determine the number of cores on your machine")
    
  }
  
  alleles <- LocusControl$alleles[loci] %>% 
    dplyr::bind_rows(.id = "locus")
  
  ploidy <- LocusControl %>% 
    dplyr::filter(locusnames %in% loci) %>% 
    dplyr::pull(ploidy) %>% 
    purrr::set_names(loci)
  
  if(sum(ploidy == 1)){
    
    stop("One or more of the loci supplied are haploid. GENEPOP only works with diploid data (ploidy = 2)")
    
  }
  
  if(is.null(npops)){
    
    np <- length(sillyvec)
    
  }else(np <- npops)
  
  sillylist <- split(sillyvec, base::ceiling(seq_along(sillyvec)/np))
  
  nalleles <- LocusControl %>% 
    dplyr::filter(locusnames %in% loci) %>% 
    dplyr::pull(nalleles) %>% 
    purrr::set_names(loci)
  
  for(i in 1:length(sillylist)){
    
    my.sillys <- sillylist[[i]]
    
    my.gcl <- sapply(sillylist[[i]], function(silly){
      
      get(paste(silly, ".gcl", sep=""), pos = 1)
      
    }, simplify = FALSE)
    
    file <- "GENEPOP input format"
    
    file <- rbind(file, cbind(loci))
    
    cl <- parallel::makePSOCKcluster(ncores)
    
    doParallel::registerDoParallel(cl, cores = ncores)  
    
    # multicore loop pop scores
    
    `%dopar%` <- foreach::`%dopar%`
    
    scores_all <- foreach::foreach(silly = my.sillys, .packages = c("tidyverse", "tidyselect")) %dopar% {
      
      new.gcl <- my.gcl[[silly]]
      
      IDs <- paste(new.gcl$SILLY_CODE, new.gcl$FK_FISH_ID, sep = "_")  # can't use SillySource since pooled pops preserve collection SillySource
      
      if(!VialNums){ 
        
        vials <- new.gcl$SILLY_CODE
        
      } else {
        
        vials <- IDs
        
      }
      
      scores_names <- sapply(loci, function(locus) {c(locus, paste0(locus, ".1"))}) %>% 
        as.vector() 
      
      scores <- new.gcl %>% 
        dplyr::select(tidyselect::all_of(scores_names)) %>% 
        as.data.frame(stringsAsFactors = FALSE)
      
      dimnames(scores)[[1]] = IDs 
      
      if(usat) {
        
        maxchar <- max(nchar(alleles$call))
        
        pop_scores <- lapply(loci, function(loc){
          
          variables <- c(loc, paste(loc, 1, sep = "."))
          
          scores %>%
            dplyr::select(tidyselect::all_of(variables)) %>%
            tidyr::replace_na(replace = list(0, 0) %>% 
                                purrr::set_names(variables)) %>%
            dplyr::mutate(dplry::across(dplyr::everything(), ~stringr::str_pad(string = ., width = maxchar, pad = "0", side = "left"))) %>% 
            tidyr::unite(col = !!rlang::as_name(loc), tidyselect::all_of(variables), sep = "")
          
        }) %>% 
          dplyr::bind_cols() %>% 
          purrr::set_names(loci) %>% 
          dplyr::mutate(ID = vials, comma = ",") %>% 
          tidyr::unite("ID_comma", c("ID", "comma"), sep = " ") %>% 
          tidyr::unite("comb_scores", tidyselect::all_of(loci), sep = " ") %>% 
          tidyr::unite("ID_comma_scores", c("ID_comma","comb_scores"), sep = "  ") %>% 
          dplyr::pull(ID_comma_scores) 
        
        rbind("Pop", cbind(pop_scores))
        
      } else {
        
        maxchar <- max(nchar(alleles$allele)) + 1
        
        pop_scores <- lapply(loci, function(loc){
          
          print(loc)
          
          variables <- c(loc, paste(loc, 1, sep = "."))
          
          my.alleles <- alleles %>% 
            dplyr::filter(locus == loc)
          
          scores %>%
            tibble::as.tibble() %>% 
            dplyr::select(tidyselect::all_of(variables)) %>% 
            dplyr::mutate(dplyr::across(dplyr::everything(), ~factor(., levels = my.alleles$call))) %>% 
            dplyr::mutate(dplyr::across(dplyr::everything(), ~as.numeric(.))) %>% 
            dplyr::mutate(dplyr::across(dplyr::everything(), ~as.character(.))) %>%
            tidyr::replace_na(replace = list("0", "0") %>%
                                purrr::set_names(variables)) %>%
            dplyr::mutate(dplyr::across(dplyr::everything(), ~stringr::str_pad(., width = maxchar, pad = "0", side = "left"))) %>% 
            tidyr::unite(col = !!rlang::as_name(loc), tidyselect::all_of(variables), sep = "")
          
        }) %>% 
          dplyr::bind_cols() %>% 
          purrr::set_names(loci) %>% 
          dplyr::mutate(ID = vials, comma = ",") %>% 
          tidyr::unite("ID_comma", c("ID", "comma"), sep = " ") %>% 
          tidyr::unite("comb_scores", tidyselect::all_of(loci), sep = " ") %>% 
          tidyr::unite("ID_comma_scores", c("ID_comma","comb_scores"), sep = "  ") %>% 
          dplyr::pull(ID_comma_scores) 
        
        rbind("Pop", cbind(pop_scores))
        
      }
      
    }  # End multicore loop
    
    parallel::stopCluster(cl)
    
    file <- rbind(file, scores_all %>% unlist() %>% cbind())
    
    popno_maxchar <- nchar(length(sillyvec))
    
    if(!is.null(npops)){ #Generic file path if npops is supplied
      
      firstpop <- stringr::str_pad(match(head(sillylist[[i]], 1), sillyvec), width = popno_maxchar, side = "left", pad = "0")

      lastpop <- stringr::str_pad(match(tail(sillylist[[i]], 1), sillyvec), width = popno_maxchar, side = "left", pad = "0")

      if(npops ==1){
        
        my.path <- paste0(path, "/Pop", firstpop, ".gen.txt")}else{
          
          my.path <- paste0(path, "/Pop", firstpop, "to", lastpop, ".gen.txt")
          
        }
      
    }else{my.path = path} #Else, keep original file path
    
    utils::write.table(x = file, file = my.path, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  print(Sys.time() - start_time)
  
  return(NULL) 
  
}