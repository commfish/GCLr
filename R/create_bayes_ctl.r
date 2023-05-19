create_bayes_ctl <- function(sillyvec, loci, mixvec, baseline_name, nreps = 40000, nchains, groupvec, priorvec, initmat, dir, seeds = matrix(sample(seq(10000), 3*nchains), nrow = 3), thin = c(1, 1, 1), mixfortran, basefortran, switches = "F T F T F T F"){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function creates a BAYES control (.ctl) file for each mixture in mixvec and each chain in nchains.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   sillyvec - character vector of population sillys used to create the baseline file
  #   loci - character vector of the loci used to produce the baseline and mixture files
  #   mixvec - a character vector of mixture ".gcl" objects you want to produce control files for. 
  #   dir - character vector of where to save the ".ctl" files
  #   baseline_name - character string of the baseline file name without the .bse extension
  #   nreps - numeric; the number of MCMC iterations you want to analyze the mixtures for.
  #   nchains - the number of MCMC chains to analyzed the mixtures for.
  #   groupvec - a numberic vector of group numbers corresponding to each silly in sillyvec
  #   priorvec - a numberic vector of priors corresponding to each silly in sillyvec
  #   initmat - a numeric matrix of intitial starting values for each chain; nrow(initmat) == length(sillyvec) & ncol(initmat) == nchains
  #   dir - the directory file path where you want the chain file to be saved
  #   seeds - a matrix of random seeds containing 3 seeds per chain; nrow(seeds) == 3 & ncol(seeds) == nchains
  #   thin - thinning intervals for MCMC sample of 1) stock proportions, 2) baseline allele or type relative frequencies, and 3) stock assignments of each mixture individual
  #   mixfortran - the fortran format of the mixture files; this is returned by create_bayes_mix
  #   basefortran - the fortran format of the baseline file; this is returned by create_bayes_base
  #   switches - a character string of logical switches; default is "F T F T F T F":
  #               1. baseline file printed, 
  #               2. mixture file printed, 
  #               3. MCMC sample of baseline relative frequencies printed, 
  #               4. baseline file contains counts, not relative frequencies, 
  #               5. stock assignment distribution for each individual in the mixture sample printed (see Stock Assignment Distribution of Mixture Individuals, pg. 12)4, 
  #               6. stock-group or regional estimates printed (see Group Estimates, pg. 11), and
  #               7. individual stock estimates are suppressed (not printed).  Options are turned “on” with “T” for true or “off” with “F” for false. 
  #             
  #  See bayes manual for addtional details:  V:\Software\BAYES\MANUAL.DOC
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Writes out bayes control (.ctl) files
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # load("V:/Analysis/2_Central/Chinook/Susitna River/Susitna_Chinook_baseline_2020/Susitna_Chinook_baseline_2020.Rdata")
  # loki2r(sillyvec = c("KSUSC18FW", "KSUSCN18", "KSUSC19FW", "KSUSCN19"), username = "awbarclay", password = password)
  # combine_loci(sillyvec = c("KSUSC18FW", "KSUSCN18", "KSUSC19FW", "KSUSCN19"), markerset = c("Ots_MHC1", "Ots_MHC2"))
  # pool_collections(collections = c("KSUSC18FW", "KSUSCN18"), newname = "Susitna2018")
  # pool_collections(collections = c("KSUSC19FW", "KSUSCN19"), newname = "Susitna2019")
  # loci <- c(loci82, "Ots_MHC1.Ots_MHC2")
  # sillyvec <- Final_Pops$silly
  # combine_loci(sillyvec = sillyvec, markerset = c("Ots_MHC1", "Ots_MHC2"))
  # base_fortran <- create_bayes_base(sillyvec = sillyvec, loci = loci, dir = "bayes/baseline", baseline_name = "SusitnaChinook31pops82loci", ncores = 8)
  # mix_fortran <- create_bayes_mix(mixvec = c("Susitna2018", "Susitna2019"), loci = loci, dir = "bayes/mixture", ncores = 8)
  # inits <- random_inits(groupvec = groupvec3, groupweights = rep(1/max(groupvec3), max(groupvec3)), nchains = 5)
  # priorvec <- create_prior(groupvec = groupvec3, groupweights = rep(1/max(groupvec3), max(groupvec3)))
  # create_bayes_control(sillyvec = sillyvec, loci = loci, mixvec = c("Susitna2018", "Susitna2019"), baseline_name = "SusitnaChinook31pops82loci", nreps = 40000, nchains = 5, groupvec = groupvec3, priorvec = priorvec, initmat = inits, dir = "bayes/control", seeds = matrix(sample(seq(10000), 3*5), nrow = 3),
  #                                  thin = c(1, 1, 1), mixfortran = mix_fortran, basefortran = base_fortran, switches = "F T F T F T F")
  # 
  # 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  if(!all(loci %in% LocusControl$locusnames)){
    
    stop(paste0("\n'", setdiff(loci, LocusControl$locusnames), "' from argument 'loci' not found in 'LocusControl' object!!!"))
    
  }
  
  sapply(mixvec, function(mix){
    
    chains <- paste("Chain", 1:nchains, sep = "")
    
    dirs <- paste(dir, "/", mix, "Chain", 1:nchains, ".ctl", sep = "") %>% setNames(chains)

    priorvec <- cbind(substring(format(round(priorvec, 6), nsmall = 6), first = 2))
    
    initmat <- cbind(substring(format(round(initmat, 6), nsmall = 6), first = 2))
    
    if(nchains > 1){
      
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
 
    files <- lapply(chains, function(chain){rbind(files[[chain]], paste0(mix, chain, "FRQ.FRQ"))})%>% setNames(chains)
    
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
      
      rbind(files[[chain]], cbind(sapply(1:nloci, function(d){paste(sprintf("%3s", d), sprintf("%2s", nalleles[loci[d]]), sprintf("%2s", ifelse(ploidy[loci[d]]==2, "T", "F")), "", loci[d], collapse = "")})))
      
      }) %>% setNames(chains)
    
    files <- lapply(chains, function(chain){
      
      rbind(files[[chain]], cbind(sapply(1:nsillys, function(i){
        
        paste(sprintf("%3s", i), sprintf("%2s", groupvec[i]), sprintf("%7s", priorvec[i]), "", format(ifelse(nchar(sillyvec[i])<18, sillyvec[i], substr(sillyvec[i], start = 1, stop = 18)), width = 18), sprintf("%7s", initmat[i, chain]), collapse = "")
        
        })))
      
      }) %>% setNames(chains)

    empty <- sapply(chains, function(chain){write.table(files[[chain]], dirs[chain], quote = FALSE, row.names = FALSE, col.names = FALSE)})
 
    }) # End mix loop   
  
}