#' @title Calculate Failure Rates by Silly, Locus, Plate, and Project
#'
#' @description
#' This function calculates failure rates by SILLY, locus, plate, and project.
#' It should be run after the "remove_na_indv" step but before any other steps
#' in the quality control (QC) process, such as removing duplicate fish,
#' alternate fish, or missing fish. Note that this function does not connect to LOKI;
#' it simply calculates failure rates (0's / total fish run) from the silly.gcl objects.
#'
#' @param sillyvec Character vector of SILLY codes in the project.
#'
#' @return A list of failure rates by SILLY, locus, plate, and project.
#'
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom plotly ggplotly
#'
#' @aliases failure_rate FailureRate.GCL.R
#'
#' @details
#'   The [GCLr::failure_rate()] function calculates failure rates by SILLY, locus, plate, and project. It performs the following steps:
#'   - Pools all collections into one master silly using the [GCLr::pool_collections()] function
#'   - Creates a tibble of Dose 1 scores and attributes from the master silly.
#'   - Calculates failure rates by SILLY, locus, plate, and project using group_by and summarise functions from the `dplyr` package.
#'   - Generates plots of failure rates by SILLY and locus, as well as by plate and locus, using `ggplot2` and `plotly` packages.
#'   - Returns a list of failure rates and plots.
#'
#'   Note: This function relies on the [GCLr::pool_collections()] function from the GCL package.
#'
#' @seealso
#'   \code{\link{remove_na_indv}}: Function to remove NA individuals
#'   \code{\link{pool_collections}}: Function to pool collections
#'
#' @family QC Functions
#'
#' @examples
#' failure_rate(sillyvec = c("SILLY1", "SILLY2", "SILLY3"))
#' 
#' @export

failure_rate <- function(sillyvec) {
  # Check for the use of the old function name
  if (match.call()[[1]] %in% c("FailureRate.GCL.R")) {
    warning("The function name 'FailureRate.GCL.R' is deprecated. Please use 'failure_rate' instead.")
  }
  
  
  # Pool all collections in to one master silly
  pool_collections(collections = sillyvec, loci = loci, newname = "master")
  
  # Tibble of Dose 1 scores and attributes
  master.tbl <- dplyr::bind_cols(as_tibble(master.gcl$scores[, , "Dose1"]), as_tibble(master.gcl$attributes[, c("SILLY_CODE", "PLATE_ID", "SillySource")])) %>% 
    tidyr::gather(locus, genotype, -SILLY_CODE, -PLATE_ID, -SillySource) %>% 
    dplyr::rename(silly = SILLY_CODE, plate = PLATE_ID, silly_source = SillySource)
  
  rm(master.gcl, pos = 1)
  
  # Failure rate by silly
  fail_silly <- master.tbl %>% 
    dplyr::group_by(silly) %>% 
    dplyr::summarise(fail = sum(genotype == "0", na.rm = TRUE) / n()) %>% 
    dplyr::arrange(dplyr::desc(fail))
  
  # Failure rate by locus
  fail_locus <- master.tbl %>% 
    dplyr::group_by(locus) %>% 
    dplyr::summarise(fail = sum(genotype == "0", na.rm = TRUE) / n()) %>% 
    dplyr::arrange(dplyr::desc(fail))
  
  # Failure rate by plate
  fail_plate <- master.tbl %>% 
    dplyr::group_by(plate) %>% 
    dplyr::summarise(fail = sum(genotype == "0", na.rm = TRUE) / n()) %>% 
    dplyr::arrange(dplyr::desc(fail))
  
  # Failure rate overall
  fail_overall <- master.tbl %>% 
    dplyr::mutate(project = project) %>% 
    dplyr::group_by(project) %>% 
    dplyr::summarise(fail = sum(genotype == "0", na.rm = TRUE) / n())
  
  # Plot failure rate by silly and locus
  fail_silly_plot <- plotly::ggplotly(
    ggplot(
      master.tbl %>%
        dplyr::group_by(silly, locus) %>%
        dplyr::summarise(p_fail = sum(genotype == "0") / n()),
      aes(x = silly, y = locus)
    ) +
      geom_tile(aes(fill = p_fail)) +
      scale_fill_gradientn(
        colours = colorRampPalette(colors = c("white", "black"))(101),
        values = seq(0.00, 1.00, by = 0.01),
        na.value = "red",
        limit = c(0, 1)
      ) +
      theme(axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )) +
      ggtitle("Failure Rate by Silly and Locus")
  )
  
  # Plot failure rate by plate and locus
  fail_plate_plot <- plotly::ggplotly(
    ggplot(
      master.tbl %>%
        dplyr::group_by(plate, locus) %>%
        dplyr::summarise(p_fail = sum(genotype == "0") / n()),
      aes(x = plate, y = locus)
    ) +
      geom_tile(aes(fill = p_fail)) +
      scale_fill_gradientn(
        colours = colorRampPalette(colors = c("white", "black"))(101),
        values = seq(0.00, 1.00, by = 0.01),
        na.value = "red",
        limit = c(0, 1)
      ) +
      theme(axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )) +
      ggtitle("Failure Rate by Plate and Locus")
  )
  
  failure_rate <- list(silly_failure_rate = fail_silly, 
                       locus_failure_rate = fail_locus, 
                       plate_failure_rate = fail_plate, 
                       overall_failure_rate = fail_overall, 
                       plot_silly_failure_rate = fail_silly_plot, 
                       plot_plate_failure_rate = fail_plate_plot)
  
  return(failure_rate)
}
