#' @title Plot Fisher's Exact Test Results
#'
#' @description 
#' This function plots the results from [GCLr::fishers_test()]. It displays a histogram of locus-specific p-values and 
#' the overall p-value across loci, faceted by each test.
#'
#' @param pooling_test Raw output from [GCLr::fishers_test()].
#' 
#' @returns A plot faceted by each test. Each facet shows:
#'     \itemize{
#'       \item a histogram of locus-specific p-values
#'       \item the overall p-value across loci printed in red
#'       \item a horizontal red line showing the expected number of loci within each p-value bin by chance
#'     }
#'
#' @details
#' This function is used to visualize the results from [GCLr::fishers_test()], which is used to test whether temporally or geographically similar 
#' baseline collections have significant differences in allele frequencies across loci. If there are no significant differences between collections 
#' (overall p-value across loci > 0.01), then we typically pool collections into populations to improve estimates of allele frequencies. However, 
#' it can also be helpful to see the distribution of locus-specific p-values to determine if there are just a few loci driving the overall p-value, or
#' if there are several loci that have low p-values.
#'
#' @examples
#' sillys <- GCLr::base2gcl(GCLr::ex_baseline, unpool = TRUE)
#' loci <- GCLr::ex_LocusControl$locusnames[-c(10,12,13,32,33)]
#' freq <- GCLr::calc_freq_pop(sillyvec = sillys, loci = loci, LocusCtl = GCLr::ex_LocusControl)
#' test.results <- GCLr::fishers_test(freq = freq, loci = loci, tests = list(c("SRAIL97","SJOHN97"), c("SMOOK93","SMOOK94"), c("SPTAR92","SPTAR93"),c("STERN92","STERN93")), LocusCtl = GCLr::ex_LocusControl)
#' GCLr::plot_fishers_test(pooling_test = test.results)
#'
#' @references While testing for homogeneity of allele frequencies is not the same as testing for conformance to Hardy-Weinberg expectations,
#' the idea for the plots created by this function comes from:  
#'   - Waples, R.S., 2015. Testing for Hardy–Weinberg proportions: have we lost the plot?. Journal of heredity, 106(1), pp.1-19. <https://academic.oup.com/jhered/article/106/1/1/2961870>
#'
#' @export
plot_fishers_test <- function(pooling_test) {
  
  # just setting variables for use in plotting - set ncol to something reasonable in facet_wrap
  if (length(pooling_test$test_sillys) > 4) {
    ncols <- round(length(pooling_test$test_sillys) / 3, 0)
  } else {
    ncols <- NULL
  }
  
  pooling_test %>%
    tidyr::unnest(bylocus) %>% {
      ggplot2::ggplot(., ggplot2::aes(x = pval)) +
        ggplot2::geom_histogram(binwidth = 0.05) +
        ggplot2::geom_hline(
          yintercept = (dplyr::select(.data = ., locus) %>% dplyr::n_distinct()) / 20,
          colour = "red"
        ) + 
        ggplot2::geom_text(
          mapping = ggplot2::aes(x = 0.5, y = 15, label = overall),
          colour = "red",
          size = 6
        ) +
        ggplot2::facet_wrap( ~ test_sillys, ncol = ncols) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "locus-specific p-values",
                      y = "number of loci")
    }
}