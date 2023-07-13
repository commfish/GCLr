#' @title QC script
#' @description Perform quality control for lab projects.
#' 
#' @note It is important to run this script in order, or some functions will not provide accurate results. This is by design, as this script only hits LOKI once to save time.
#' @note This is not a function, but a step wise operation that the user must run through to manually generate applicable files.
#' 
#' @param dirqc The directory where the QC files are located.
#' @param species The species for the project (e.g., "pink").
#' @param project The project name (e.g., "P014").
#' @param project_type The type of project ("Biomark", "uSat", or "GT-seq").
#' @param username The username for LOKI access.
#' @param password The password for LOKI access.
#' @param qcSummaryfile The file path for the QC summary file; This is automatically generated using `project` (don't change unless you know what you're doing).
#' @param conflict_rate The conflict rate threshold at which dupcheck between sillys occurs (default is 0.10).
#'
#' @returns Returns summary tables and saves them in the QC summary file.
#' \itemize{
#'   \item \code{text}: DupCheck Results (if applicable)
#'   \item \code{figures}: 
#'    - Levelplot of genotyping success rate by silly and locus for ProjectSillys
#'    - Levelplot of genotyping success rate by silly and locus for qcSillys
#'    - Histogram of overall qc individual conflict rate
#'    - Individual histograms of duplicate rate for conflict individuals
#'   \item \code{summary excel}: 
#'    - Summary by Silly
#'    - Conflicts by Silly
#'    - DupCheck Results (if applicable)
#'    - Conflicts by Locus
#'    - Conflicts by PlateID
#'    - Failure Rate by Silly
#'    - Failure Rate by Locus
#'    - Failure Rate by PlateID
#'    - Overall Failure Rate
#'    - Original Project Sample Size by Locus
#'    - Duplicates within Silly
#' }
#' 
#' @examples
#' \dontrun{
#' # Example usage
#' dirqc <- "V:/Lab/Genotyping/SNP Projects/Pink/Project P014 AHRP Parentage GTSeq Part 2/qc/"
#' species <- "pink"
#' project <- "P014"
#' project_type <- "Biomark"
#' username <- ""
#' password <- ""
#' qcSummaryfile <- paste("Project", project,"qc Summary R Script.xlsx")
#' conflict_rate <- 0.10
#' P014_qc(dirqc, species, project, project_type, username, password, qcSummaryfile, conflict_rate)
#' }
#'
#' @export
qc <- function(dirqc, species, project, project_type, username, password, qcSummaryfile = paste("Project", project,"qc Summary R Script.xlsx"), conflict_rate = 0.10) {
  #stop it from running automatically
  if (FALSE) {
   
    # Title: ProjectTitle
    # Date: TodaysDate
    # Name: YourName
    

    #~~~  GO! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    setwd(dirqc)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Read in Project Genotypes ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    GCLr::read_proj_geno(project_name = project, username = username, password = .password)
    
    rm(.password)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Failure Rate ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    (failure_rate <- GCLr::failure_rate(sillyvec = ProjectSillys))
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Read in qc Genotypes ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    qcfiles <- list.files(path = "Genotype Data Files", pattern = ".csv", full.names = TRUE, recursive = FALSE)
    
    if (project_type == "Biomark") {
      GCLr::read_biomark_qc.r(qccsvFilepaths = qcfiles)
    } else {
      if (project_type == "uSat") {
        GCLr::read_usat_qc(qccsvFilepaths = qcfiles)
      } else {
        if (project_type == "GT-seq") {
          GCLr::read_gtseq_qc(qccsvFilepaths = qcfiles)
        } else {
          message("project_type must be either 'Biomark', 'uSat', or' GT-seq' to read in qc genotypes")
        }  # GT-seq
      }  # uSat
    }  # Biomark
    
    qcColSize <- sapply(paste(qcSillys, ".gcl", sep = ''), function(x) get(x)$n)
    
    qcColSizeAll <- setNames(rep(0, length(ProjectSillys)),paste0(ProjectSillys, "qc.gcl"))
    
    qcColSizeAll[paste0(qcSillys, ".gcl")] <- qcColSize[paste0(qcSillys, ".gcl")]
    
    qcColSizeAll
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Read in Conflict Report ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    qcConcordanceReportfile <- list.files(path = "Conflict Reports", pattern = "ConcordanceReport", full.names = TRUE)
    
    GCLr::combine_conflicts(files = qcConcordanceReportfile)
    
    # Old conflict report has "0" for mitochondrial conflicts, new has " " for mitochondrial conflicts, we will refer to them as "Homo-Homo".
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Conflict summaries
    
    conflicts_by_plate <- combined_conflicts %>% 
      dplyr::group_by(plate_id, concordance_type) %>%
      dplyr::summarise(n = n()) %>% 
      tidyr::spread(concordance_type, n, fill = 0, drop = FALSE) %>% 
      dplyr::mutate(Conflict = sum(`Het-Het`, `Het-Homo`, `Homo-Het`, `Homo-Homo`)) %>% 
      dplyr::ungroup()
    
    conflicts_by_silly <- combined_conflicts %>% 
      dplyr::group_by(silly, concordance_type) %>% 
      dplyr::summarise(n = n()) %>% 
      tidyr::spread(concordance_type, n, fill = 0, drop = FALSE) %>% 
      dplyr::mutate(Conflict = sum(`Het-Het`, `Het-Homo`, `Homo-Het`, `Homo-Homo`)) %>% 
      dplyr::ungroup()
    
    conflicts_by_locus <- combined_conflicts %>% 
      dplyr::group_by(locus, concordance_type) %>% 
      dplyr::summarise(n = n()) %>% 
      tidyr::spread(concordance_type, n, fill = 0, drop = FALSE) %>% 
      dplyr::mutate(Conflict = sum(`Het-Het`, `Het-Homo`, `Homo-Het`, `Homo-Homo`)) %>% 
      dplyr::ungroup()
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Sample Size by Locus for Project Genotypes ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    OriginalProjectSampleSizebyLocus <- GCLr::sampn_by_locus(sillyvec = ProjectSillys, loci = loci)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Sample Size by Locus for qc Genotypes ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    OriginalqcSampleSizebyLocus <- GCLr::sampn_by_locus(sillyvec = qcSillys, loci = loci)
    
    OriginalqcPercentbyLocus <- apply(OriginalqcSampleSizebyLocus, 1, function(row) {row / max(row)} )
    
    rerunsqc <- which(apply(OriginalqcPercentbyLocus, 2, min) < 0.8)
    
    
    plotly::ggplotly(
      ggplot2::ggplot(
        OriginalqcPercentbyLocus %>%
          tidyr::as_tibble(rownames = "locus") %>%
          tidyr::gather(silly, percent, -locus),
        aes(x = silly, y = locus)
      ) +
        ggplot2::geom_tile(aes(fill = percent), color = "white") +
        ggplot2::scale_fill_gradient(low = "black", high = "white", limits = c(0, 1)) +
        ggplot2::ylab("Locus") +
        ggplot2::xlab("SILLY") +
        ggplot2::ggtitle("Percent genotyped by silly/locus") +
        ggplot2::labs(fill = "Percent\ngenotyped") +
        ggplot2::theme(axis.text.x = element_text(angle = 90))
    )
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### QA of Project Genotypes ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ProjectSillys_SampleSizes <- matrix(data = NA, nrow = length(ProjectSillys), ncol = 5, dimnames = list(ProjectSillys, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))
    
    ProjectSillys_SampleSizes[, "Genotyped"] <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)
    
    if (species %in% c("chum", "sockeye")) {
      
      Alternate <- GCLr::find_alt_species(sillyvec = ProjectSillys, species = species) %>% 
        dplyr::as_tibble()
      
      nAltBySilly <- sapply(ProjectSillys, function(silly) {
        AlternateSpeciesReport <- Alternate[grep(pattern = silly, x = rownames(Alternate)), ]
        sum(AlternateSpeciesReport$Alternate > 0.5 & AlternateSpeciesReport$Failure > 0.5)
      })
      # remove_alt_species(AlternateSpeciesReport = Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)  # Do not remove fish, just note how many per silly. Still want to catch them in conflicts later.
      
    } else {
      
      Alternate = tibble::tibble(x = "Not applicable")
      
    }
    
    ColSizePostAlternate <- ProjectSillys_SampleSizes[, "Genotyped"]
    if (exists(x = "nAltBySilly")) {ColSizePostAlternate <- ColSizePostAlternate - nAltBySilly}
    # ColSizePostAlternate <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)
    
    ProjectSillys_SampleSizes[, "Alternate"] <- ProjectSillys_SampleSizes[, "Genotyped"] - ColSizePostAlternate 
    
    MissLoci <- GCLr::remove_ind_miss_loci(sillyvec = ProjectSillys, proportion = 0.8)
    
    ColSizePostMissLoci <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n) - ProjectSillys_SampleSizes[, "Alternate"]
    
    ProjectSillys_SampleSizes[, "Missing"] <-  ColSizePostAlternate - ColSizePostMissLoci
    
    DuplicateCheck95MinProportion <- GCLr::dupcheck_within_silly(sillyvec = ProjectSillys, loci = loci, quantile = NULL, minproportion = 0.95)
    
    DuplicateCheckReportSummary <- sapply(ProjectSillys, function(x) DuplicateCheck95MinProportion[[x]]$report, simplify = FALSE)
    
    nDupsBySilly <- sapply(DuplicateCheckReportSummary, function(silly) {ifelse(is.character(silly), 0, nrow(as.matrix(silly)))})
    # RemovedDups <- remove_dups(DuplicateCheck95MinProportion)  # Do not remove fish, just note how many per silly. Still want to catch them in conflicts later.
    
    sapply(DuplicateCheckReportSummary[nDupsBySilly >= 1], function(silly) {if (1 %in% abs(as.numeric(levels(silly$ID1)) - as.numeric(levels(silly$ID2)))) {"Sequential IDs found as duplicates, check 'DuplicateCheckReportSummary' for duplicated rows"} else {"Duplicates exist, but IDs do not appear sequential"} } )
    
    DuplicateCheckReportSummary[nDupsBySilly >= 1]  # Show within silly duplicates
    
    ColSizePostDuplicate <- ColSizePostMissLoci - nDupsBySilly
    # ColSizePostDuplicate <- sapply(paste(ProjectSillys, ".gcl", sep = ''), function(x) get(x)$n)
    
    ProjectSillys_SampleSizes[, "Duplicate"] <- ColSizePostMissLoci - ColSizePostDuplicate
    
    ProjectSillys_SampleSizes[, "Final"] <- ColSizePostDuplicate
    
    ProjectSillys_SampleSizes
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### QA of qc Genotypes ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    MissLociqc <- GCLr::remove_ind_miss_loci(sillyvec = qcSillys, proportion = 0.8)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Perform Duplicate Check on High Conflict Individuals ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Filter for conflicts, determine conflict rate
    conflicts <- combined_conflicts %>% 
      dplyr::filter(concordance_type %in% c("Het-Het", "Het-Homo", "Homo-Het", "Homo-Homo")) %>% 
      dplyr::count(silly_source) %>% 
      dplyr::mutate(p = n / length(loci))
    
    # Histogram of conflict rate
    conflicts %>% 
      ggplot2::ggplot(aes(x = p)) +
      ggplot2::geom_bar() +
      ggplot2::xlim(0, 1) +
      ggplot2::geom_vline(xintercept = conflict_rate, colour = "red", lwd = 1.5) +
      ggplot2::xlab("Conflict rate") +
      ggplot2::ylab("Frequency") +
      ggplot2::ggtitle("qc individual conflict rate")
    
    # Filter for conflicts > conflict_rate
    conflicts_investigate <- conflicts %>% 
      dplyr::filter(p > conflict_rate)
    
    # Duplicate check if necessary
    if (nrow(conflicts_investigate) == 0) {
      
      message(paste0("No individuals have > ", conflict_rate * 100, "% loci with conflicts between project and qc."))
      
      dup_check_results <- tibble::tibble(x = "Not applicable")
      
    } else {
      
      message(paste0("The following individuals have > ", conflict_rate * 100, "% loci with conflicts between project and qc:\n"), paste(conflicts_investigate$silly_source, conflicts_investigate$n, "conflicts", collapse = "\n"))
      
      # Loop through individuals to see if missing loci
      conflict_indv <- NULL
      
      for (silly_ind in conflicts_investigate$silly_source) {
        
        silly <- stringr::str_split(string = silly_ind, pattern = "_", simplify = TRUE)[, 1]
        
        ind <- stringr::str_split(string = silly_ind, pattern = "_", simplify = TRUE)[, 2]
        
        # qc fish lost in QA?
        if (ind %in% MissLociqc[[paste0(silly, "qc")]]) {
          message(paste0("\n", silly, "qc_", ind, " does not have at least 80% loci genotyped, not running DupCheck for this individual."))
        }  # if qc fish removed due to missing genotypes
        
        # Project fish lost in QA
        if (ind %in% MissLoci[[silly]]) {
          message(paste0("\n", silly, "_", ind, " does not have at least 80% loci genotyped, not running DupCheck for this individual."))
        }  # if project fish removed due to missing genotypes
        
        conflict_indv <- c(conflict_indv, paste(silly, ind, sep = "_")[!(ind %in% MissLociqc[[paste0(silly, "qc")]] | ind %in% MissLoci[[silly]]) ])  # Confirm qc fish and Project fish were not removed
        
      }  # silly_ind
      
      # If no more, stop
      if (is.null(conflict_indv) | length(conflict_indv) == 0) {
        
        message("\nNo remaining high conflict individuals.")
        
        dup_check_results <- tibble::tibble(x = "Not applicable")
        
      } else {
        
        conflicts_investigate <- conflicts_investigate %>% 
          dplyr::filter(silly_source %in% conflict_indv)
        
        message("\nRunning dupcheck_among_sillys on these high conflict individuals, as they have at least 80% loci genotyped for Project and qc extractions.")
        message(paste(conflicts_investigate$silly_source, conflicts_investigate$n, "conflicts", collapse = "\n"))
        
        conflict_silly <- unique(stringr::str_split(string = conflicts_investigate$silly_source, pattern = "_", simplify = TRUE)[, 1])
        
        KeySillyIDs <- setNames(
          lapply(conflict_silly, function(silly) {
            sapply(grep(pattern = silly, x = conflict_indv, value = TRUE), function(ind) {
              stringr::str_split(string = ind, pattern = "_", simplify = TRUE)[, 2]
            }, USE.NAMES = FALSE) 
          }),
          paste0(conflict_silly, "qc"))
        
        DupCheckResults <- sapply(conflict_silly, function(silly) {
          GCLr::dupcheck_among_sillys(KeySillys = paste0(silly, "qc"), 
                                KeySillyIDs = KeySillyIDs[paste0(silly, "qc")], 
                                BetweenSillys = ProjectSillys, 
                                loci = loci, 
                                threshold = 0.9)
        }, simplify = FALSE)  # FALSE
        
        dup_check_results <- dplyr::bind_rows(DupCheckResults, .id = "silly") %>% 
          tibble::as_tibble()
        
      }  # conflict_ind, post missing individuals
      
    }  # else, conflicts_to_investigate
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Create Summary Tables ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    summary_table_1 <- dplyr::bind_cols(tibble(Silly = ProjectSillys), as_tibble(ProjectSillys_SampleSizes)) %>% 
      dplyr::left_join(failure_rate$silly_failure_rate, by = c("Silly" = "silly")) %>% 
      dplyr::rename("Failure Rate" = fail) %>% 
      dplyr::mutate("Total qc Fish" = qcColSizeAll)
    
    qc_silly_genotypes <- tibble::tibble(silly = factor(ProjectSillys),
                                         qc_genotypes = sapply(ProjectSillys, function(silly) {
                                           qc_silly = paste0(silly, "qc.gcl")
                                           ifelse(qc_silly %in% names(qcColSizeAll), qcColSizeAll[qc_silly] * length(loci), 0)
                                         } ))
    
    summary_table_2 <- conflicts_by_silly %>% 
      tidyr::gather(type, number, -silly) %>%  # make tall
      dplyr::left_join(qc_silly_genotypes, by = "silly") %>%  # join number of qc genotypes by silly
      dplyr::mutate(rate = number / qc_genotypes) %>%  # conflict numbers to rates
      tidyr::gather(variable, value, -silly, -qc_genotypes, -type) %>%  # make tall
      tidyr::unite(temp, type, variable) %>%  # unite conflict type with both number and rate
      tidyr::spread(temp, value) %>%  # make wide
      dplyr::rename(Silly = silly, 
                    "Total qc Genotypes" = qc_genotypes, 
                    "Total Discrepancies" = Conflict_number, 
                    "Discrepancy Rate" = Conflict_rate,
                    "DB Zeros" = `DB Zero_number`,
                    "DB Zero Rate" = `DB Zero_rate`,
                    "qc Zeros" = `File Zero_number`,
                    "qc Zero Rate" = `File Zero_rate`,
                    "Total Het-Het" = `Het-Het_number`,
                    "Het-Het Rate" = `Het-Het_rate`,
                    "Total Het-Homo" = `Het-Homo_number`,
                    "Het-Homo Rate" = `Het-Homo_rate`,
                    "Total Homo-Het" = `Homo-Het_number`,
                    "Homo-Het Rate" = `Homo-Het_rate`,
                    "Total Homo-Homo" = `Homo-Homo_number`,
                    "Homo-Homo Rate" = `Homo-Homo_rate`)
    
    summary_table_3 <- conflicts_by_locus %>% 
      tidyr::gather(type, number, -locus) %>%  # make tall
      dplyr::mutate(qc_genotypes = sum(qcColSizeAll)) %>%  # join number of qc genotypes by locus
      dplyr::mutate(rate = number / qc_genotypes) %>%  # conflict numbers to rates
      tidyr::gather(variable, value, -locus, -qc_genotypes, -type) %>%  # make tall
      tidyr::unite(temp, type, variable) %>%  # unite conflict type with both number and rate
      tidyr::spread(temp, value) %>%  # make wide
      dplyr::rename(Locus = locus,
                    "Total qc Genotypes" = qc_genotypes, 
                    "Total Discrepancies" = Conflict_number, 
                    "Discrepancy Rate" = Conflict_rate,
                    "DB Zeros" = `DB Zero_number`,
                    "DB Zero Rate" = `DB Zero_rate`,
                    "qc Zeros" = `File Zero_number`,
                    "qc Zero Rate" = `File Zero_rate`,
                    "Total Het-Het" = `Het-Het_number`,
                    "Het-Het Rate" = `Het-Het_rate`,
                    "Total Het-Homo" = `Het-Homo_number`,
                    "Het-Homo Rate" = `Het-Homo_rate`,
                    "Total Homo-Het" = `Homo-Het_number`,
                    "Homo-Het Rate" = `Homo-Het_rate`,
                    "Total Homo-Homo" = `Homo-Homo_number`,
                    "Homo-Homo Rate" = `Homo-Homo_rate`)
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Append Summary Tables to qcSummaryfile.xlsx ####
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Create list object for all output to simple summary.xlsx file
    summary_lst <- suppressWarnings(list("Summary by Silly" = summary_table_1,
                                         "Conflicts by Silly" = summary_table_2,
                                         "Conflicts by Locus" = summary_table_3,
                                         "Conflicts by PlateID" = conflicts_by_plate,
                                         "qc Duplicate Check" = dup_check_results,
                                         "Failure Rate by Silly" = failure_rate$silly_failure_rate,
                                         "Failure Rate by Locus" = failure_rate$locus_failure_rate,
                                         "Failure Rate by Plate" = failure_rate$plate_failure_rate,
                                         "Overall Failure Rate" = failure_rate$overall_failure_rate,
                                         "Project Sample Size by Locus" = OriginalProjectSampleSizebyLocus %>% tibble::rownames_to_column("silly") %>% tibble::as_tibble(),
                                         "Duplicate Check in Project" = dplyr::bind_rows(DuplicateCheckReportSummary[!DuplicateCheckReportSummary == "No Duplicates"], .id = "silly") %>% tibble::as_tibble(),
                                         "Alternate Species" = Alternate))
    
    # Write out a "Simple" file, can't update normal Summary File by inserting new tabs
    if (file.exists(qcSummaryfile)) {
      stop(paste0("qc Summary file: '", qcSummaryfile ,"' already exists, change the file name so you don't overwrite it, hoser!!!"))
    } else {
      writexl::write_xlsx(x = summary_lst, path = qcSummaryfile, col_names = TRUE)
    }
    
    #~~~  STOP!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
#stop from running automatically
  }
}

