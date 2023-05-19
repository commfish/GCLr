plot_SampSizeByLocus <- function(SampSizeByLocus){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function takes the output from sampn_by_locus() and creates a plotly heatmap of the proportion of fish with scores for each locus and silly
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   SampSizeByLocus - a tibble of sample sizes by locus produced by sampn_by_locus(). 
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  
  #  A plotly heatmap of the proportion of fish with scores for each locus and silly
  #
  # Examples~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  create_locuscontrol(markersuite = "Sockeye2011_96SNPs", username = "awbarclay", password = password)
  #  sillyvec = c("SMCDO03", "SNEVA13")
  #  password = "************"
  #  loki2r(sillyvec = sillyvec, username = "awbarclay", password = password)
  #  remove_ind_miss_loci(sillyvec = sillyvec)
  #
  #  Output <- sampn_by_locus(sillyvec = sillyvec, loci = LocusControl$locusnames)
  #  plot_SampSizeByLocus(SampSizeByLocus = Output)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  sillyvec <- SampSizeByLocus$silly
  
  silly_n <- silly_n(sillyvec)
  
  plotly::ggplotly(
    
    SampSizeByLocus %>% 
      tidyr::pivot_longer(-silly, names_to = "locus", values_to = "count") %>% 
      dplyr::full_join(silly_n, by = "silly") %>% 
      dplyr::mutate(proportion = count/n) %>% 
      ggplot2::ggplot(aes(x = silly, y = locus, fill = proportion))+
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low = "white", high = "blue4", limits = c(0, 1)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
    
  )
  
}