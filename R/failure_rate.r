failure_rate <- function(sillyvec) {
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  This function calculates failure rates by silly, loci, plate, and project
#
#  It is IMPORTANT to run this function after "remove_na_indv" but before any
#  other steps in the qc (i.e. removing duplicate fish, alternate fish, 
#  missing fish, etc.). This function does NOT connect to LOKI. It merely
#  calculates failure rates (0's / total fish run) from the silly.gcl objects.
#
#  Argument(s):  
#  sillyvec <- character vector of SILLYs in the project
#
#  Output:
#  List of failure rate by SILLY, locus, and project
#
#  Written by Kyle Shedd 10/16/15  
#  Updated for tidyverse by Kyle Shedd 9/5/18
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  while(!require(tidyverse)){ install.packages("tidyverse") }
  
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
      geom_tile (aes(fill = p_fail)) +
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
      geom_tile (aes(fill = p_fail)) +
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