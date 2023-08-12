#' Convert GENEPOP File to COLONY Input File
#'
#' Create COLONY genotype input file from a GENEPOP file.
#'
#' @param path_genepop The full file path of the GENEPOP file with either .gen or .txt extension.
#' @param path_colony The full file path of the COLONY file including .txt extension.
#' 
#' @details 
#' This function is designed to convert a GENEPOP file to a COLONY genotype file. 
#' If you are doing parentage analysis, you'll need to run this separately for parents and offspring.
#' If you are using adult sex data, you'll run the parents separately for Mum and Dad.
#' 
#' @return A .txt file of genotypes for COLONY. The first column is the ID, followed by a 2 column genotypes (2 columns per locus).
#' 
#' @examples
#' GCLr::gcl2genepop(sillyvec = "unuk_parr", loci = LocusControl$locusnames, path = "../data/genotypes/combined_postQA/unuk_parr.gen", VialNums = TRUE, usat = FALSE, ncores = 4)
#' GCLr::genepop2colony(path_genepop, path_colony)
#' @export
genepop2colony <- function(path_genepop, path_colony) {
  rawdat <- scan(path_genepop, what = "", sep = "\n")
  
  len <- length(rawdat)
  
  popind <- sapply(rawdat, function(lin) {
    strsplit(lin, " ")[[1]][1] == "Pop"
    
  }) %>% as.vector()
  
  npops <- sum(popind)
  
  ORD <- order(popind, decreasing = TRUE)[1:npops]
  
  loci <- sapply(rawdat[2:(ORD[1] - 1)], function(str) {
    strsplit(str, " ")[[1]][1]
    
  })
  
  nloc <- length(loci)
  
  dat0 <-
    rawdat[-(1:(nloc + 1))][rawdat[-(1:(nloc + 1))] != "Pop"] %>%
    tibble::as_tibble() %>%
    tidyr::separate(
      col = value,
      sep = " ,  ",
      into = c("SillySource", "geno")
    ) %>%
    tidyr::separate(col = geno, sep = " ", into = loci)
  
  
  numeric_geno <- dat0 %>%
    tidyr::pivot_longer(cols = -1,
                        names_to = "locus",
                        values_to = "genotype") %>%
    tidyr::separate(
      col = "genotype",
      into = c("allele_1", "allele_2"),
      sep = 3,
      remove = TRUE
    ) %>%
    dplyr::mutate(allele_1 = as.numeric(allele_1),
                  allele_2 = as.numeric(allele_2))
  
  allele_1 <- numeric_geno %>%
    dplyr::select(SillySource, locus, allele_1) %>%
    tidyr::pivot_wider(names_from = locus, values_from = allele_1)
  
  allele_2 <- numeric_geno %>%
    dplyr::select(SillySource, locus, allele_2) %>%
    tidyr::pivot_wider(names_from = locus, values_from = allele_2)
  
  dat_geno <-
    dplyr::left_join(
      x = allele_1,
      y = allele_2,
      by = "SillySource",
      suffix = c("", ".1")
    ) %>%
    dplyr::select(SillySource, tidyselect::all_of(paste0(
      rep(LocusControl$locusnames, each = 2), c("", ".1")
    )))
  
  readr::write_delim(
    x = dat_geno,
    file = path_colony,
    delim = "\t",
    col_names = FALSE
  )
}