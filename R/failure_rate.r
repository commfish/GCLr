#' @title Calculate Failure Rates by Silly, Locus, Plate, and Overall
#'
#' @description
#' This function calculates failure rates by SILLY, locus, plate, and overall.
#' The main purpose for this function is for use in the [GCLr::qc()] markdown; however, it can also be use in other analyses (i.e., baseline, mixture) to calculate failure rates.
#' Note that this function does not connect to LOKI;
#' it simply calculates failure rates (0's / total fish run) from the silly.gcl objects.
#'
#' @param sillyvec Character vector of SILLY codes you want to include in the failure rate calculation.
#' @param loci A character vector of locus names to include in the failure rate calculation (default = LocusControl$locusnames)
#'
#' @returns A list of failure rates by SILLY, locus, plate, and overall.
#'
#' @details
#'   The [GCLr::failure_rate()] function calculates failure rates by SILLY, locus, plate, and project. It performs the following steps:
#'   - Pools all collections into one master silly using the [GCLr::pool_collections()] function
#'   - Creates a tibble of Dose 1 scores and attributes from the master silly.
#'   - Calculates failure rates by SILLY, locus, plate, and overall using [dplyr::group_by()] and [dplyr::summarise()].
#'   - Generates plots of failure rates by SILLY and locus, as well as by plate and locus, using `ggplot2` and `plotly` packages.
#'   - Returns a list of failure rates and plots.
#'
#'   Note: This function relies on the [GCLr::pool_collections()] function from the GCL package.
#'
#' @seealso
#'   \code{\link{pool_collections}}: Function to pool collections
#'
#' @family QC Functions
#'
#' @examples
#'   
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' loci <- GCLr::ex_baseline[,-c(1:5)] %>%
#'   names() %>%
#'   gsub(pattern = "*\\.1", x = ., replacement = "") %>%
#'   unique()
#' 
#' GCLr::failure_rate(sillyvec = sillyvec, loci = loci, LocusCtl = GCLr::ex_LocusControl)
#' 
#' @export
failure_rate <- function(sillyvec, loci = LocusControl$locusnames, LocusCtl = LocusControl) {
  # Pool all collections in to one master silly
  GCLr::pool_collections(collections = sillyvec,
                         loci = loci,
                         newname = "master",
                         LocusCtl = LocusCtl)
  
  # Tibble of Dose 1 scores and attributes
  master.tbl <- master.gcl %>%
    dplyr::mutate(silly = stringr::str_remove(SillySource, "_.*")) %>% # grab silly from silly_source
    dplyr::select(silly,
                  plate = PLATE_ID,
                  silly_source = SillySource,
                  tidyselect::all_of(loci)) %>%
    tidyr::pivot_longer(
      cols = tidyselect::starts_with(loci),
      names_to = "locus",
      values_to = "genotype"
    ) %>%
    dplyr::arrange(locus, silly_source)
  
  rm(master.gcl, pos = 1)
  
  # Add in a fake plate ID if none exists. This is needed to plot failure rate by plate and locus
  if (any(is.na(master.tbl$plate))) {
    
    message("Some or all of the .gcl objects supplies do not contain plate IDs. Assigning a fake plate ID ('00000') where they are missing.")
    
    master.tbl <- master.tbl %>% 
      dplyr::mutate(plate = dplyr::case_when(is.na(plate)~"00000",
                                      TRUE~plate))

  }
   
  # Failure rate by silly
  fail_silly <- master.tbl %>%
    dplyr::group_by(silly) %>%
    dplyr::summarise(fail = sum(is.na(genotype), na.rm = FALSE) / dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(fail))
  
  # Failure rate by locus
  fail_locus <- master.tbl %>%
    dplyr::group_by(locus) %>%
    dplyr::summarise(fail = sum(is.na(genotype), na.rm = FALSE) / dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(fail))
  
  # Failure rate by plate
  fail_plate <- master.tbl %>%
    dplyr::group_by(plate) %>%
    dplyr::summarise(fail = sum(is.na(genotype), na.rm = FALSE) / dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(fail))
  
  # Failure rate overall
  fail_overall <- master.tbl %>%
    dplyr::summarise(fail = sum(is.na(genotype), na.rm = FALSE) / dplyr::n())
  
  # Plot failure rate by silly and locus
  fail_silly_plot <- plotly::ggplotly(
    ggplot2::ggplot(
      master.tbl %>%
        dplyr::group_by(silly, locus) %>%
        dplyr::summarise(
          p_fail = sum(is.na(genotype), na.rm = FALSE) / dplyr::n(),
          .groups = "drop"
        ),
      ggplot2::aes(x = silly, y = locus)
    ) +
      ggplot2::geom_tile(ggplot2::aes(fill = p_fail)) +
      ggplot2::scale_fill_gradientn(
        colours = grDevices::colorRampPalette(colors = c("white", "black"))(101),
        values = seq(0.00, 1.00, by = 0.01),
        na.value = "red",
        limit = c(0, 1)
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )) +
      ggplot2::ggtitle("Failure Rate by Silly and Locus")
  )
  
  # Plot failure rate by plate and locus
  fail_plate_plot <- plotly::ggplotly(
    ggplot2::ggplot(
      master.tbl %>%
        dplyr::group_by(plate, locus) %>%
        dplyr::summarise(
          p_fail = sum(is.na(genotype), na.rm = FALSE) / dplyr::n(),
          .groups = "drop"
        ),
      ggplot2::aes(x = plate, y = locus)
    ) +
      ggplot2::geom_tile(ggplot2::aes(fill = p_fail)) +
      ggplot2::scale_fill_gradientn(
        colours = grDevices::colorRampPalette(colors = c("white", "black"))(101),
        values = seq(0.00, 1.00, by = 0.01),
        na.value = "red",
        limit = c(0, 1)
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )) +
      ggplot2::ggtitle("Failure Rate by Plate and Locus")
  )
  
  failure_rate <- list(
      silly_failure_rate = fail_silly,
      locus_failure_rate = fail_locus,
      plate_failure_rate = fail_plate,
      overall_failure_rate = fail_overall,
      plot_silly_failure_rate = fail_silly_plot,
      plot_plate_failure_rate = fail_plate_plot
    )
  
  return(failure_rate)
}