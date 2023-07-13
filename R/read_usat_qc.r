#' @title Read uSat QC Genotypes
#'
#' @description This function reads in uSat QC genotypes from .csv files and creates .gcl objects in the global environment. It verifies that the QC sillys and loci are found in the project and assigns the `qcSillys` vector to the global environment. It's called on by [qc()].
#'
#' @param qccsvFilepaths A character vector with relative paths for QC genotype CSV files.
#'
#' @details The function reads in the QC genotype data from the CSV files, renames and formats columns, and creates scores and counts arrays for each silly. The resulting .gcl objects are assigned to the global environment with the names "sillyqc.gcl" (e.g., "silly1qc.gcl", "silly2qc.gcl").
#'
#' @returns Returns a few silly objects to the global environment 
#'  `qc.gcl` list objects
#'  `qcSillys`; a character vector of qc sillys
#'   
#' @export
read_usat_qc <- function(qccsvFilepaths) {
  # Read in .csv files
  qc_genotypes <- suppressMessages(
    suppressWarnings(
      dplyr::bind_rows(
        lapply(qccsvFilepaths, function(fle) {readr::read_csv(file = fle)[, 1:4]} )
      )  
    )  
  ) 
  
  # Rename columns, split silly_source
  qc_genotypes <- qc_genotypes %>% 
    dplyr::rename(silly_source = "Sample Name", locus = Marker, allele_1 = "Allele 1", allele_2 = "Allele 2") %>% 
    tidyr::separate(col = silly_source, into = c("silly", "fish_id"), sep = "_", remove = FALSE) %>% 
    dplyr::mutate(fish_id = as.character(as.numeric(fish_id))) %>% 
    tidyr::unite("silly_source", c(silly, fish_id), sep = "_", remove = FALSE)
  
  # Verify that all qc silly are in the project
  ProjectSillysqc <- unique(qc_genotypes$silly)
  if (!all(ProjectSillysqc %in% ProjectSillys)) { stop(paste0(ProjectSillysqc[!ProjectSillysqc %in% ProjectSillys], " not found in ProjectSillys.")) }
  
  # Verify that all qc loci are in project loci
  lociqc <- sort(unique(qc_genotypes$locus))
  if (!all(lociqc %in% loci)) { stop(paste0(lociqc[!lociqc %in% loci], " not found in LocusControl.")) }
  
  # attributes table names
  attNames <- colnames(get(paste0(ProjectSillys[1], ".gcl"))$attributes)
  
  # Loop over silly to create .gcl objects
  for (x in ProjectSillysqc) {
    
    # subset genotypes by silly
    x_genotypes <- qc_genotypes %>% 
      dplyr::filter(silly == x) %>% 
      dplyr::mutate(locus = factor(locus, levels = loci))
    
    # create scores array
    scores <- Reduce(bbind, lapply(c("allele_1", "allele_2"), function(al){ tapply(dplyr::pull(x_genotypes, al), list(x_genotypes$fish_id, x_genotypes$locus), c)[,, drop = FALSE] }))
    dimnames(scores)[[3]] <- paste0("Dose", 1:2)
    
    # create counts array
    counts <- array(NA, c(nrow(scores), ncol(scores), max(nalleles)), list(rownames(scores), colnames(scores), paste0("Allele", seq(max(nalleles)))))
    for (locus in loci) {
      for (al in seq(nalleles[locus])) {
        for (id in x_genotypes$fish_id) {
          counts[id, locus, al] <- sum(scores[id, locus, seq(ploidy[locus])] == alleles[[locus]][al])
        }  # id
      }  # al           
    }  # locus
    
    # create attributes data.frame
    attributes <- data.frame(matrix(NA, nrow = nrow(counts), ncol = length(attNames), dimnames = list(rownames(counts), attNames)))
    attributes$FK_FISH_ID <- rownames(counts)
    attributes$SillySource <- paste0(x, "qc_", rownames(counts))
    
    # assign to global environment
    assign(x = paste0(x, "qc.gcl"), value = list(counts = counts, scores = scores, n = nrow(scores), attributes = attributes), pos = 1)
    
  }  # x (silly)
  
  assign(x = "qcSillys", value = paste0(ProjectSillysqc, "qc"), pos = 1)
}
