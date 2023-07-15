#' @title Plot Genepop HWE Output
#'
#' @description 
#' This function plots the results from [GCLr::read_genepop_hwe()]. It displays a histogram of locus- or population-specific p-values and 
#' the overall p-value across loci/populations, faceted by each locus/population.
#' 
#' @param GenepopHWE_report Raw output from [GCLr::read_genepop_hwe()].
#' @param sillyvec An optional character vector of silly codes without the ".gcl" extension used to filter results when `plot_type = "silly"` (default = `NULL`).
#' If `NULL`, `sillyvec` will be inherited from `GenepopHWE_report`. This argument does not filter results when `plot_type = "locus"`.
#' @param plot_type A Character vector of length 1 specifying whether to plot p-values by "silly" or "locus" (default = "silly").
#' If "locus", only loci with an overall p-value across populations < 0.1 will be shown.
#'
#' @returns A plot faceted by each population/locus Each facet shows:
#'     \itemize{
#'       \item a histogram of locus- or population-specific p-values
#'       \item the overall p-value across loci/populations printed in red
#'       \item a horizontal red line showing the expected number of loci/populations within each p-value bin by chance
#'     }
#'
#' @details
#' This function is used to visualize the results from [GCLr::read_genepop_hwe()], which is used to test whether collections/populations have significant 
#' deviations from Hardy-Weinberg expectations (HWE; overall p-value < 0.05 with Bonferroni correction for the number of collections in the baseline).
#' When using this function to test loci for conformance to HWE, we recommend removing loci when the overall p-value < 0.01. However, in both cases
#' it can also be helpful to see the distribution of locus- or population-specific p-values to determine if there are just a few loci/populations driving the 
#' overall p-value, or if there are several loci/populations that have low p-values. When populations or loci do not conform to HWE, users should investigate
#' patterns of Fis to determine the root cause.
#' 
#' @examples
#' \dontrun{
#' genepop_hwe <- read_genepop_hwe(file = "~/R/test.txt.P", sillyvec = sillyvec)
#' plot_genepop_hwe(GenepopHWE_report = genepop_hwe, plot_type = "silly")
#' plot_genepop_hwe(GenepopHWE_report = genepop_hwe, plot_type = "locus")
#' }
#' 
#' @references The idea for the plots created by this function comes from:  
#'   - Waples, R.S., 2015. Testing for Hardy–Weinberg proportions: have we lost the plot?. Journal of heredity, 106(1), pp.1-19. <https://academic.oup.com/jhered/article/106/1/1/2961870>
#'   
#' @export
plot_genepop_hwe <- function(GenepopHWE_report, sillyvec = NULL, plot_type = c("silly", "locus")[1]) {

  summary_sillys <- setdiff(dimnames(GenepopHWE_report$SummaryPValues)[[2]], "Overall Pops")
  
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
  } else {
    sillyvec <- summary_sillys
  }
  
  if (plot_type == "silly") {
    
    # just setting variables for use in plotting - set ncol to something reasonable in facet_wrap
    if (length(sillyvec) > 4) {
      ncols <- round(length(sillyvec) / 3, 0)
    } else {
      ncols <- NULL
    }
    
    # First convert the Genepop HWE report into a usable table and filter for sillyvec.
    HWEpval <- GenepopHWE_report$SummaryPValues %>%
      dplyr::as_tibble(rownames = "locus") %>%
      tidyr::pivot_longer(-locus, names_to = "silly", values_to = "pval")  %>%
      dplyr::filter(silly %in% sillyvec)
    
    # Now plot the pvals
    HWEpval %>% {
      ggplot2::ggplot(dplyr::filter(.data = ., locus != "Overall Loci"), ggplot2::aes(x = pval)) +
        ggplot2::geom_histogram(binwidth = 0.05) +
        ggplot2::geom_hline(
          yintercept = (
            dplyr::filter(.data = ., locus != "Overall Loci") %>% dplyr::select(locus) %>% dplyr::n_distinct()
          ) / 20,
          colour = "red"
        ) +
        ggplot2::geom_text(
          data = dplyr::filter(.data = ., locus == "Overall Loci"),
          mapping = ggplot2::aes(x = 0.5, y = 15, label = pval),
          colour = "red",
          size = 6
        ) +
        ggplot2::facet_wrap( ~ silly, ncol = ncols) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "locus-specific p-values",
                      y = "number of loci")
    }
    
  } else if (plot_type == "locus") {

    # First convert the Genepop HWE report into a usable table and filter for loci.
    HWEpval <- 
      GenepopHWE_report$SummaryPValues %>%
      dplyr::as_tibble(rownames = "locus") %>%
      tidyr::pivot_longer(-locus, names_to = "silly", values_to = "pval")
    
    low_loci <-
      HWEpval %>% 
      filter(silly == "Overall Pops" &
               pval < 0.1, locus != "Overall Loci") %>% 
      pull(locus)
    
    # just setting variables for use in plotting - set ncol to something reasonable in facet_wrap
    if (length(low_loci) > 4) {
      ncols <- round(length(low_loci) / 3, 0)
    } else {
      ncols <- NULL
    }
    
    # Now plot the pvals by batches according to ncols:
    HWEpval %>% {
      ggplot2::ggplot(dplyr::filter(.data = ., silly != "Overall Pops" &
                                      locus %in% low_loci),
                      ggplot2::aes(x = pval)) +
        ggplot2::geom_histogram(binwidth = 0.05) +
        ggplot2::geom_hline(
          yintercept = (dplyr::select(.data = ., silly) %>% dplyr::n_distinct()) / 20,
          colour = "red"
        ) +
        ggplot2::geom_text(
          data = dplyr::filter(.data = ., silly == "Overall Pops" &
                                 locus %in% low_loci),
          mapping = ggplot2::aes(x = 0.5, y = 15, label = pval),
          colour = "red",
          size = 6
        ) +
        ggplot2::facet_wrap(~ locus, ncol = ncols) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "pop-specific p-values",
                      y = "number of pops")
    }
    
  } else {
    stop("Did you specify plotting by `silly` or by `locus`?")
  }
}