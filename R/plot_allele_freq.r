plot_allele_freq <- function(freq, file, sillyvec = NULL, groupvec = NULL, loci = NULL, group_col = NULL, popnames = NULL, xlab.cex = 6, ylab.cex = 12, xtick.int = 1){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function uses the output from calc_freq_pop() and produces allele freqency bubble plots.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   freq - a frequency tibble produced by calc_freq_pop()
  #
  #   file - the file path including the file name with .pdf extention where the plots will be saved
  #   
  #   sillyvec - a vector of silly codes without the ".gcl" extention (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART10")). 
  #
  #   groupvec - numeric vector indicating the group affiliation of each pop in sillyvec
  #
  #   group_col - vector of colors corresponding to groupvec
  #
  #   loci - vector of locus names; if set to NULL all loci in the ".gcl" obejects will be used.
  #
  #   popnames - a vector of new population names (e.g. popnames <- c("KQUART","KCUP","KPINT"); Default is NULL, which assigns numbers seq_along(sillyvec) to popnames. 
  #              Note: Setting popnames for a large number of pops may cause the x-axis tick labels to become unreadable
  #
  #   xlab.cex  - text size on x axis
  # 
  #   ylab.cex - text size y axis
  #
  #   xtick.int - x-axis tick interval. This argument only works when popnames = NULL. Default is set to 1, which will produce ticks and tick labels for all pops.  
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  A pdf file containing allele fequency plots for each locus in loci
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  load("V:\\Analysis\\5_Coastwide\\Chum\\NPen2WA_Chum_baseline\\NPen2WA_Chum_baseline.Rdata")
  #  Freq <- calc_freq_pop(sillyvec = sillyvec227, loci = loci91)
  #   
  #  plot_allele_freq(freq = Freq, file = "./test_freq_plots.pdf", sillyvec = sillyvec227, groupvec = groupvec19, loci = loci91, group_col = grcol, popnames = NULL, xlab.cex = 6, ylab.cex = 12, xtick.int = 10)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  if(is.null(loci)) {
    loci = freq$locus %>% unique()
  }
  
  if (is.null(groupvec)) {
    
    groupvec = seq(length(freq$silly %>% unique()))
    
  }
  
  if (is.null(group_col)) {
    
    group_col = rainbow(max(groupvec))
    
  }
  
  if(is.null(sillyvec)) {
    
    sillyvec = freq$silly %>% unique()
    
  }
  
  if(!length(sillyvec)==length(groupvec)) {
    
    stop("Make sure sillyvec and groupvec are the same length!!")
    
  }
  
  # set popnames to numeric if NULL
  if (is.null(popnames)) {
    
    popnames = seq_along(sillyvec)
    
  } else{
    
    warning("Setting popnames for a large number of pops may cause the x-axis tick labels to become unreadable. 
           If this occurs, try setting popnames = NULL to label them with population numbers and increase xtick.int(x tick interval) to reduce the number of ticks and tick labels on the x-axis.")
    
    }
  
  # Making it esier to read numbers on the x-axis of plots
  if (is.numeric(popnames)) {
    
    angle = 45
    
  } else {
    
    angle = 90
    
  }
  
  freq$proportion[freq$proportion == 0] <- NA 
  
  freq_df <- freq %>% 
    dplyr::left_join(dplyr::tibble(silly = sillyvec, groupvec), 
                     by = "silly") %>% 
    dplyr::mutate(groupvec = factor(groupvec))
  
  if(is.numeric(popnames)){
    
    brks <- c(1, seq(0, length(popnames), by = xtick.int))
    
    breaks <- sillyvec[brks]
    
    labels <- factor(popnames)[brks]
      
  } else{
    
    breaks <- ggplot2::waiver()
    
    labels <- popnames
    
  }
  
  pdf(file = file)
  
  lapply(loci, function(locus){
    
    plot <- freq_df %>%
      dplyr::filter(locus == locus) %>% 
      dplyr::filter(locus == !!locus) %>% 
      dplyr::mutate(silly = factor(silly, levels = sillyvec), 
                    allele = factor(allele, unique(allele))) %>% 
      ggplot2::ggplot(aes(x = silly, y = allele, color = groupvec, size = proportion*100)) +
      ggplot2::geom_point() + 
      ggplot2::scale_color_manual(values = group_col, guide = FALSE) +
      ggplot2::scale_x_discrete(labels = labels, breaks = breaks) +
      ggplot2::xlab("Pop")+
      ggplot2::theme(legend.position = "none", 
                     axis.text.x = ggplot2::element_text(angle = angle, hjust = 1, vjust = 0.5, size = xlab.cex), 
                     axis.text.y = ggplot2::element_text(size = ylab.cex)) +
      ggplot2::ggtitle(label = locus)
    
   suppressWarnings(print(plot))
    
    })
  
  dev.off()

}