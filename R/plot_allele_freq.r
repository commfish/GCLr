#' @title Plot Allele Frequency Bubble Plots
#'
#' @description 
#' This function plots the results from [GCLr::calc_freq_pop()]. It creates a ("*.pdf") of bubble plots of allele frequencies for each locus.
#' 
#' @param freq A tibble of allele frequencies generated by [GCLr::calc_freq_pop()].
#' @param file A character vector of length 1 with the full file path including ".pdf" extension where plots will be saved.
#' @param sillyvec An optional character vector of silly codes without the ".gcl" extension (default = `NULL`).
#' If `NULL`, `sillyvec` will come directly from `freq`.
#' @param groupvec An optional numeric vector indicating the group affiliation of each silly in `sillyvec`, 
#' must be the same length as `sillyvec` (default = `NULL`).
#' @param loci An optional character vector of locus names (default = `NULL`).
#' If `NULL`, `loci` will come directly from `freq`.
#' @param group_col An optional character vector of colors corresponding to each group in `groupvec` (default = `NULL`).
#' @param popnames An optional character vector of collection/population names for each silly in `sillyvec`,
#' must be the same length as `sillyvec` (default = `NULL`).
#' If `NULL`, `popnames` will be `1:length(sillyvec)`.
#' Setting `popnames` for a large number of pops may cause the x-axis tick labels to become unreadable.
#' @param xlab.cex An optional numeric vector of length 1 indicating the x-axis label font size, (default = 6).
#' @param ylab.cex An optional numeric vector of length 1 indicating the y-axis label font size, (default = 12).
#' @param xtick.int An optional numeric vector of length 1 indicating the x-axis tick interval, (default = 1).
#' This argument only works when `popnames = NULL`.
#'
#' @returns A ("*.pdf") file containing allele frequency bubble plots, 1 page per locus
#'
#' @details
#' Allele frequency bubble plots are most useful for loci with multiple alleles (e.g. uSATs or microhaplotypes). 
#' For visualizing allele frequencies for SNPs, use [GCLr::plot_freq_fis_4snps()].
#' 
#' @seealso 
#' [GCLr::calc_freq_pop()]
#' [GCLr::plot_freq_fis_4snps()]
#' 
#' @examples
#' \dontrun{
#' freqs <- GCLr::calc_freq_pop(sillyvec = sillyvec)
#' GCLr::plot_allele_freq(freq = freqs, file = "~/R/test_allele_freq.pdf")
#' }
#' 
#' @export
plot_allele_freq <-
  function(freq,
           file,
           sillyvec = NULL,
           groupvec = NULL,
           loci = NULL,
           group_col = NULL,
           popnames = NULL,
           xlab.cex = 6,
           ylab.cex = 12,
           xtick.int = 1) {

  if (is.null(loci)) {
    loci = freq$locus %>% unique()
  }
  
  if (is.null(groupvec)) {
    
    groupvec = seq(length(freq$silly %>% unique()))
    
  }
  
  if (is.null(group_col)) {
    
    group_col = rainbow(dplyr::n_distinct(groupvec))
    
  }
  
  if (is.null(sillyvec)) {
    
    sillyvec = freq$silly %>% unique()
    
  }
  
  if (!length(sillyvec) == length(groupvec)) {
    
    stop("Make sure sillyvec and groupvec are the same length!!")
    
  }
  
  # set popnames to numeric if NULL
  if (is.null(popnames)) {
    
    popnames = seq_along(sillyvec)
    
  } else {
    
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
  
  if (is.numeric(popnames)){
    
    brks <- c(1, seq(0, length(popnames), by = xtick.int))
    
    breaks <- sillyvec[brks]
    
    labels <- factor(popnames)[brks]
      
  } else {
    
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
      ggplot2::ggplot(ggplot2::aes(x = silly, y = allele, color = groupvec, size = proportion * 100)) +
      ggplot2::geom_point() + 
      ggplot2::scale_color_manual(values = group_col, guide = FALSE) +
      ggplot2::scale_x_discrete(labels = labels, breaks = breaks) +
      ggplot2::xlab("Pop") +
      ggplot2::theme_bw() + 
      ggplot2::theme(legend.position = "none", 
                     axis.text.x = ggplot2::element_text(angle = angle, hjust = 1, vjust = 0.5, size = xlab.cex), 
                     axis.text.y = ggplot2::element_text(size = ylab.cex)) +
      ggplot2::ggtitle(label = locus)
    
   suppressWarnings(print(plot))
    
    })
  
  dev.off()

}