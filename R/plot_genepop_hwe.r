plot_genepop_hwe <- function(GenepopHWE_report, sillyvec = NULL, plot_type = "silly") {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function plots results from read_genepop_hwe(). You provide the report, a list of sillys, and loci and it will provide a visual of
  #   p-values faceted by silly.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   GenepopHWE_report - raw output from read_genepop_hwe()
  #
  #   sillyvec <- vector of sillys you're interested in
  #
  #   plot_type <- do you want to show plots by silly or by locus? Options are"silly" or "locus"; default = "silly"
  #                 if by locus, only the low-pvalue loci will be shown.
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Produces a plot of the results from read_genepop_hwe. Specifically, displays p-value, overall p-value, by silly or by locus.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  #  attach("V:\\Analysis\\5_Coastwide\\Chum\\NPen2WA_Chum_baseline\\NPen2WA_Chum_baseline.Rdata")
  #  sillyvec234 <- sillyvec234
  #  detach()
  #  
  #  HWEreport <- read_genepop_hwe("V:/Analysis/5_Coastwide/Chum/NPen2WA_Chum_baseline/GENEPOP/NAKPen2WA_234pops_93loci.P", sillyvec = sillyvec234)
  #
  #  plot_genepop_hwe(GenepopHWE_report = HWEreport, sillyvec = sillyvec234[1:10])
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (!require("pacman"))
    install.packages("pacman")
  library(pacman)
  pacman::p_load(tidyverse)  # Install packages, if not in library and then load them
  
  summary_sillys <- dimnames(GenepopHWE_report$SummaryPValues)[[2]]
  
  # Sillyvec is optional in read_genepop_hwe.
  # Check to make sure all sillys have results if sillyvec is supplied.
  if (!is.null(sillyvec)) {
    n_sillys <- length(sillyvec)
    
    n_matches <-  sum(summary_sillys %in% sillyvec)
    
    if (n_sillys > n_matches) {
      stop(
        "GeneopHWE_report does not contain results for all sillys in sillyvec. Was sillyvec supplied for read_genepop_hwe?"
      )
      
    }
    
  } else{
    sillyvec <- summary_sillys
    
  }
  
  # just setting variables for use in plotting - set ncol to something reasonable in facet_wrap
  if (length(sillyvec) > 4) {
    ncols <- round(length(sillyvec) / 3, 0)
    
  } else {
    ncols <- NULL
    
  }
  
  if (plot_type %in% "silly") {
    # First convert the Genepop HWE report into a usable table and filter for sillyvec.
    HWEpval <- GenepopHWE_report$SummaryPValues %>%
      dplyr::as_tibble(rownames = "locus") %>%
      tidyr::pivot_longer(-locus, names_to = "silly", values_to = "pval")  %>%
      dplyr::filter(silly %in% sillyvec)
    
    # Now plot the pvals
    ggplot2::ggplot(data = HWEpval %>% dplyr::filter(locus != "Overall Loci"), aes(x = pval)) +
      ggplot2::geom_histogram(binwidth = 0.05) +
      ggplot2::geom_hline(
        yintercept = (HWEpval %>% dplyr::select(locus) %>% dplyr::n_distinct()) / 20,
        colour = "red") +
      ggplot2::geom_text(
        data = HWEpval %>% dplyr::filter(locus == "Overall Loci"),
        mapping = aes(x = 0.5, y = 15, label = pval),
        colour = "red",
        size = 6) +
      ggplot2::facet_wrap( ~ silly, ncol = ncols) +
      ggplot2::theme_bw()
    
  } else if (plot_type %in% "locus") {
    # First convert the Genepop HWE report into a usable table and filter for loci.
    HWEpval <- GenepopHWE_report$SummaryPValues %>%
      dplyr::as_tibble(rownames = "locus") %>%
      tidyr::pivot_longer(-locus, names_to = "silly", values_to = "pval")
    
    low_loci <-
      HWEpval %>% filter(silly == "Overall Pops" &
                           pval < 0.1, locus != "Overall Loci") %>% pull(locus)
    
    # Now plot the pvals by batches according to ncols:
    ggplot2::ggplot(data = HWEpval %>% 
                      dplyr::filter(silly != "Overall Pops" & locus %in% low_loci),
                    aes(x = pval)) +
      ggplot2::geom_histogram(binwidth = 0.05) +
      ggplot2::geom_hline(
        yintercept = (HWEpval %>% dplyr::select(silly) %>% dplyr::n_distinct()) / 20,
        colour = "red") +
      ggplot2::geom_text(
        data = HWEpval %>% dplyr::filter(silly == "Overall Pops" &
                                           locus %in% low_loci),
        mapping = aes(x = 0.5, y = 15, label = pval),
        colour = "red",
        size = 6) +
      ggplot2::facet_wrap( ~ locus, ncol = round(length(low_loci) / 10, 0)) +
      ggplot2::theme_bw()
    
  } else {
    stop("Did you specify plotting by silly or by locus?")
    }
}