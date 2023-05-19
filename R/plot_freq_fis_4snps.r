plot_freq_fis_4snps <- function(sillyvec, loci, groupvec, alpha = 0.05, groupcol = NULL, file = NULL, group.pch = 19, point.size = 1, pval.cex = 3, pval.digits = 2, line.width = .5, ncores = 4){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function will create a pdf file with plots of allele frequency and Fis 
  # for each locus on a separate page. HWE p-values are printed if less than
  # alpha.
  # 
  # Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  #  sillyvec - a vector of sillys, without ".gcl" extension that you want the 
  #   freqencies plotted for
  #
  #  loci - a vector of loci that you want included in the frequency plots, 
  #   **Note: Only use loci with ploidy==2. If loci with ploidy == 1 are supplied, 
  #           they will not be included in the plots and a message will appear.**
  #
  #  groupvec - a vector of group numbers corresponding to each silly in sillyvec
  #
  #  alpha - critical HWE p-value
  #
  #  groupcol - optional character vector of colors corresponding to each group of length max(groupvec), 
  #             can be either color numbers (numeric vector) or color names (character vector).
  #             If left NULL (default) a rainbow of colors will be used.
  #     
  #  file - the full file path, with .pdf extension where the file will be written
  #         If no file is supplied, the default is to write the file "FreqPlot.pdf" to 
  #         the current working directory.
  #
  #  group.pch - a vector of pch values corresponding to each group with length = 
  #              max(groupvec). If only one pch is supplied it will be recycled for each group; default pch is 19
  #
  #  point.size - the size of points on the plot, used by ggplot2::geom_point()
  #
  #  pval.cex - character expansion (i.e. size of text) for p-values
  #
  #  line.width - the width of the line connecting the points, used by ggplot2::geom_line()
  #
  #  ncores - the number of cores to use in a foreach %dopar% loop. If the number of core exceeds the number on your device, then ncores defaults to parallel::detectCores() 
  # 
  # Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  A pdf file with allele frequency and Fis plots for
  #  each locus and indicates the HWE pvalue if less than alpha.
  #
  #  Returns a tibble with the following variables:
  #   1) population
  #   2) locus
  #   3) pval
  #   4) fis
  #   5) allele_freq
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  load("V:\\Analysis\\5_Coastwide\\Chum\\NPen2WA_Chum_baseline\\NPen2WA_Chum_baseline.Rdata")
  #   
  #  plot_freq_fis_4snps(sillyvec = sillyvec227, loci = loci91, groupvec = groupvec19, alpha = 0.05, groupcol = grcol, file = "./FreqFisPlot.pdf", ncores = 8)
  #    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  start.time <- Sys.time()

  if(is.null(file)){

    file <- paste0(getwd(),"/FreqFisPlot.pdf")

  }
  
  
  if(ncores > parallel::detectCores()) {ncores = parallel::detectCores()}

  ploidy <- LocusControl$ploidy[loci]

  nalleles <- LocusControl$nalleles[loci]
  
  if(length(setdiff(loci, loci[ploidy==2 & nalleles==2]))>=1){
    
    message(paste0("The following loci were excluded from the plots becuase they have a ploidy of 1: ", paste0(setdiff(loci, loci[ploidy==2 & nalleles==2]), collapse = ", ")))
    
    }
  
  loci <- loci[ploidy==2 & nalleles==2]

  get_counts <- function(my.gcl, loci){
    
    lapply(loci, function(locus){
      
      my.alleles <- LocusControl$alleles[[locus]] %>% 
        pull(call)
      
      my.gcl %>% 
        dplyr::select(contains(locus)) %>% 
        dplyr::mutate_all(.funs = factor, levels = my.alleles) %>% 
        dplyr::mutate_all(.funs = as.numeric) %>% 
        dplyr::mutate_all(.funs = as.character) %>%
        dplyr::mutate(counts = rowSums(.=="1")) %>% 
        dplyr::pull(counts)
    
    }) %>% 
      dplyr::bind_cols() %>% 
      purrr::set_names(loci)

  }
  
  gcls <- lapply(sillyvec, function(silly){get(paste0(silly,".gcl"), pos = 1)}) %>% 
    purrr::set_names(sillyvec)

  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  #Start parallel loop
  HWE_df <- foreach::foreach(silly = sillyvec, .export = "LocusControl", .packages = c("tidyverse", "HardyWeinberg")) %dopar% {

    my.gcl <- gcls[[silly]]

    genecounts0 <- 3 - get_counts(my.gcl, loci) 

    genecounts <- t(apply(genecounts0, 2, tabulate, nbins = 3))

    HWEpval <- suppressWarnings(HardyWeinberg::HWExactMat(genecounts, verbose = FALSE))$pvalvec

    FixBOOL <- apply(genecounts, 1, function(gcounts){sum(gcounts==0)>1})

    HWEfis <- suppressWarnings(base::apply(genecounts[!FixBOOL, ], 1, HardyWeinberg::HWf)) %>% 
      tibble::as_tibble() %>% 
      dplyr::rename(fis = value) %>% 
      dplyr::mutate(locus = dimnames(genecounts[!FixBOOL, ])[[1]])
   
    tibble::tibble(pop_no = match(silly, sillyvec) %>% na.omit() %>% as.numeric(), population = silly, locus = loci, pval = HWEpval) %>% 
      dplyr::left_join(HWEfis, by = "locus") %>% 
      dplyr::mutate(fis = tidyr::replace_na(fis, 0))
    
  } %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(population = factor(population, levels = sillyvec))#silly
  
  parallel::stopCluster(cl) #End parallel loop

  G <- max(groupvec)

  if(is.null(groupcol)){
    
    groupcol <- rainbow(G)
    
  }
    
  if(is.numeric(groupcol)){
    
     groupcol <- colors()[groupcol]
     
  }  

  PopCol <- groupcol[groupvec] 
  
  Freq <- suppressMessages(calc_freq_pop(sillyvec = sillyvec, loci = loci, ncores = ncores))
    
  q <- Freq %>% 
    dplyr::filter(allele_no == 1) %>% 
    dplyr::pull(proportion)
  
  HWE_df$q <- q
  
  if(length(group.pch)==1){

    group.pch <- rep(group.pch, G)

  }
  
  PopPch <- group.pch[groupvec]
  
  fisylim0 <- max(abs(HWE_df$fis), na.rm = TRUE)

  fisylim <- c(-fisylim0, fisylim0)

  pdf(file, width = 11, height = 8.5, family = "Helvetica", pointsize = 12)

  for(locus in loci){
  
    my.dat <- HWE_df %>% 
      dplyr::filter(locus==!!locus)
    
    freqplot <- my.dat %>% 
      ggplot2::ggplot(aes(y = q, x = pop_no)) +
      ggplot2::geom_point(color = PopCol, shape = PopPch, size = point.size)+
      ggplot2::geom_line(color = "black", linetype = "dashed", size = line.width) +
      ggplot2::ylim(0,1)+
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank()) + 
      ggplot2::ggtitle(label = unique(locus))
    
    labels <- as.character(round(my.dat$pval , pval.digits))
    
    labels[labels=="0"] <- paste(c("<0.", rep(0, pval.digits), "5"), collapse = "")
    
    labels[my.dat$pval>alpha] <- NA
    
    lab_adjust <- fisylim0*0.1
    
    fisplot <- my.dat %>% 
      dplyr::mutate(labels = !!labels) %>% 
      ggplot2::ggplot(aes(y = fis, x = pop_no)) +
      ggplot2::geom_abline(slope = 0, size = 0.5)+
      ggplot2::geom_point(color = PopCol, shape = PopPch, size = point.size)+
      ggplot2::geom_line(color = "black", linetype = "dashed", size = line.width) +
      ggplot2::ylim(fisylim)+
      ggplot2::ylab(expression(italic(F)[IS]))+
      ggplot2::xlab("Population")+
      ggplot2::geom_text(aes(y = fis+lab_adjust*sign(fis), label = labels), color = PopCol, size = pval.cex) +
      ggplot2::theme_bw()
    
    grid::grid.newpage()
    grid::grid.draw(rbind(ggplot2::ggplotGrob(freqplot), ggplot2::ggplotGrob(fisplot), size = "last"))

  }

  dev.off()
  
  print(Sys.time()-start.time)

  return(HWE_df %>% dplyr::select(population, locus, pval, fis, allele_freq = q))
 
}

  
  