#' @name generate_haplotypes_for_loki
#' @title Generate all possible haplotypes for GT-seq amplicons
#' @description This function takes a data frame containing SNP information and 
#' generates all possible haplotype combinations for GT-seq amplicons.
#' 
#' @param input: A data frame containing the following columns:
#'   * SNP (optional)
#'   * LOCUS_NAME: The identifier/name for each locus.
#'   * SNPpos: Used for sorting, first SNP in the amplicon, then next, etc.
#'   * Allele1: The first allele at a given SNP, already alphabetical
#'   * Allele2: The second allele at a given SNP, already alphabetical
#' @param output: Path to the output CSV file. Note that it must contain the 
#' filename and extension as well (e.g., ~/../Desktop/output_haplos.csv)
#' 
#' @return A data frame containing the following columns:
#'   * LOCUS: Locus name
#'   * VALUE: Combined alleles for each haplotype
#' @examples
#' # Example Usage
#' my_haplotypes <- generate_haplotypes_for_loki(input = "path/to/input.csv", output = "path/to/output.csv")
#' @export

generate_haplotypes_for_loki <- function(input, output) {
  # Step 1: Read input data
  haps <- readr::read_csv(input, show_col_types = FALSE)
suppressMessages({
  # Step 2: Process haplotypes for each locus
  all_haps <- pbapply::pbsapply(unique(haps$LOCUS_NAME), function(locus) {
    
    # Filter and arrange data for the current locus
    x <- haps %>% 
      dplyr::filter(LOCUS_NAME == locus) %>%
      dplyr::arrange(SNPpos)  # Assuming SNPpos defines the order of alleles
    
    # Calculate parameters for haplotype generation 
    nloci <- nrow(x)
    ploidy <- 2
    
    # Generate haplotype combinations 
    hap_tib <- tibble::as_tibble(t(combn(x = rep(seq(ploidy), nloci), m = nloci)),.name_repair = "unique") %>% 
      dplyr::distinct() %>% 
      dplyr::arrange_all()  # Ensure consistent haplotype ordering 
    
    # Calculate haplotypes by combining alleles
    locus_haps <- sapply(seq(nrow(hap_tib)), function(i) {
      hap <- as.matrix(hap_tib[i, ])
      
      tmp <- sapply(seq(nloci), function(loc) {
        x[loc, hap[loc] + 3]  # Extract alleles based on haplotype 
      }) %>%
        dplyr::bind_cols()
      colnames(tmp) <- paste0("Locus", seq(nloci))
      tmp
    }, simplify = FALSE) %>%
      dplyr::bind_rows() %>%
      tidyr::unite("VALUE", sep = '') %>% 
      dplyr::add_row("VALUE" = "0") %>% 
      dplyr::mutate(LOCUS = locus) %>% 
      dplyr::select(LOCUS, VALUE)
    
    locus_haps
  }, simplify = FALSE) %>% 
    dplyr::bind_rows()
  })
  
  # Step 3: Write output CSV file
  readr::write_csv(all_haps, output)
}