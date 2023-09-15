#' Plot Sample Size by Locus
#'
#' This function takes a data frame containing information about sample sizes by locus
#' and creates a heatmap plot to visualize the proportion of counts for each locus
#' relative to the overall sample size. It uses the `plotly` and `ggplot2` packages
#' for plotting.
#'
#' @param SampSizeByLocus A data frame with two columns: "silly" and "count". The "silly"
#'   column contains categorical data, and the "count" column contains the count of samples
#'   for each category. This data frame is produced by [GCLr::sampn_by_locus()].
#'
#' @return A `plotly` interactive plot showing a heatmap of the proportion of counts for each
#'   locus relative to the overall sample size.
#'
#' @details
#' The function first calculates the proportion of counts for each locus by dividing the
#' count for each locus by the overall sample size. It then creates a long-format data frame
#' using the `pivot_longer` function from the `tidyr` package. The long-format data frame
#' includes three columns: "silly" (categorical data), "locus" (locus names), and "proportion"
#' (proportion of counts for each locus). Next, the function performs a full join between the
#' original data and the `silly_n` data, where `silly_n` is a function that computes the
#' overall sample size. The function finally uses the `ggplot2` package to create a heatmap
#' plot with "silly" on the x-axis, "locus" on the y-axis, and the "proportion" mapped to
#' the fill color.
#'
#' @examples
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' sampn_df <- GCLr::sampn_by_locus(sillyvec, loci = GCLr::ex_LocusControl$locusnames[-c(10, 12, 13, 32, 33)], LocusCtl = GCLr::ex_LocusControl)
#' 
#' GCLr::plot_sampn_by_locus(sampn_df)
#' @export
plot_sampn_by_locus <- function(SampSizeByLocus){
  
  sillyvec <- SampSizeByLocus$silly
  
  silly_n <- silly_n(sillyvec)
  
  plotly::ggplotly(
    
    SampSizeByLocus %>% 
      tidyr::pivot_longer(-silly, names_to = "locus", values_to = "count") %>% 
      dplyr::full_join(silly_n, by = "silly") %>% 
      dplyr::mutate(proportion = count/n) %>% 
      ggplot2::ggplot(ggplot2::aes(x = silly, y = locus, fill = proportion))+
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low = "white", high = "blue4", limits = c(0, 1)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
    
  )
  
}