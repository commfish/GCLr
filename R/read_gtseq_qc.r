ReadGTseqqc.GCL <- function(qccsvFilepaths) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  This function reads in GTseq qc genotypes from .csv files as .gcl objects
  #
  #  Argument(s):  
  #  qccsvFilepaths <- character vector with relative path for qc genotype .csv files
  #
  #  Output:
  #  Creates qc.gcl list objects in global environment
  #  qcSillys - a character vector of qc sillys is assigned to the global environment
  #
  #  Written by Kyle Shedd 10/15/18  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  while(!require(tidyverse)){ install.packages("tidyverse") }
  
  # Read in .csv files
  qc_genotypes <- suppressMessages(
    suppressWarnings(
      dplyr::bind_rows(
        lapply(qccsvFilepaths, function(fle) {readr::read_csv(file = fle, na = c("", "NA", "0", "0/0"))} )  # GENOTYPE == "0" is NA
      )  # bind_rows
    )  # supressWarnings
  )  # suppressMessages
  
  # Rename columns, split genotype, unite silly_source
  qc_genotypes <- qc_genotypes %>% 
    dplyr::rename(silly = SILLY_CODE, fish_id = SAMPLE_NUM, locus = LOCUS) %>% 
    dplyr::mutate(fish_id = as.character(fish_id)) %>% 
    tidyr::separate(GENOTYPE, into = c("allele_1", "allele_2"), sep = "/") %>% 
    tidyr::unite(silly_source, c(silly, fish_id), remove = FALSE) %>% 
    dplyr::select(silly_source, silly, fish_id, locus, allele_1, allele_2)
  
  # Rename columns, split silly_source
  # qc_genotypes <- qc_genotypes %>% 
  #   dplyr::rename(silly_source = "Sample Name", locus = Marker, allele_1 = "Allele 1", allele_2 = "Allele 2") %>% 
  #   tidyr::separate(col = silly_source, into = c("silly", "fish_id"), sep = "_", remove = FALSE)

  # Verify that all qc silly are in the project
  ProjectSillysqc <- unique(qc_genotypes$silly)
  if(!all(ProjectSillysqc %in% ProjectSillys)){ stop(paste0(ProjectSillysqc[! ProjectSillysqc %in% ProjectSillys], " not found in ProjectSillys.")) }
  
  # Verify that all qc loci are in project loci
  lociqc <- sort(unique(qc_genotypes$locus))
  if(!all(lociqc %in% loci)){ stop(paste0(lociqc[! lociqc %in% loci], " not found in LocusControl.")) }
  
  # attributes table names
  attNames <- colnames(get(paste0(ProjectSillys[1], ".gcl"))$attributes)
  
  # Loop over silly to create .gcl objects
  for(x in ProjectSillysqc){
    
    # subset genotypes by silly
    x_genotypes <- qc_genotypes %>% 
      dplyr::filter(silly == x) %>% 
      dplyr::mutate(locus = factor(locus, levels = loci))
    
    # create scores array
    scores <- Reduce(bbind, lapply(c("allele_1", "allele_2"), function(al){ tapply(pull(x_genotypes, al), list(x_genotypes$fish_id, x_genotypes$locus), c)[,, drop = FALSE] }))
    dimnames(scores)[[3]] <- paste0("Dose", 1:2)
    
    # create counts array
    counts <- array(data = NA, 
                    dim = c(nrow(scores), ncol(scores), max(nalleles)), 
                    dimnames = list(rownames(scores), colnames(scores), paste0("Allele", seq(max(nalleles))))
    )

    for(locus in loci){
      for(al in seq(nalleles[locus])){
        for(id in unique(x_genotypes$fish_id)){
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