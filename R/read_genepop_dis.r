#' @title Read Genepop DIS Output
#' 
#' @description
#' This function reads in output from a `genepop` Linkage Disequilibrium test ("*.DIS") file.
#'
#' @param file The full file path to the `genepop` LD output, including the ".DIS" extension.
#' @param loci An optional character vector of locus names (default = `NULL`).
#' If `NULL` the locus names will come directly from the (.DIS) file and get truncated to 8 characters.
#' If supplying `loci`, make sure it is the same `loci` used in [GCLr::gcl2genepop()].
#' @param sillyvec An optional character vector of silly codes without the ".gcl" extension (default = `NULL`).
#' If `NULL`, `Pop` names will come directly from the (".DIS") file and likely include "_fishID" extensions.
#' If supplying `sillyvec`, make sure it is the same `sillyvec` used in [GCLr::gcl2genepop()].
#'
#' @returns A non-tidy tibble with pairwise loci LD p-values per population and overall populations:
#'     \itemize{
#'       \item \code{Locus1}: locus name for 1st locus in pair
#'       \item \code{Locus2}: locus name for 2nd locus in pair
#'       \item \code{sillyvec_1}: pairwise loci LD p-values for 1st silly in `sillyvec`
#'       \item \code{sillyvec_2}: pairwise loci LD p-values in `sillyvec`
#'       \item \code{sillyvec_n}: pairwise loci LD p-values in `sillyvec`
#'       \item \code{Overall}: pairwise loci LD p-values overall sillys in `sillyvec`       
#'       
#' @details
#' Designed and tested on Genepop v4.8.3. If the LD (".DIS) file contains tests for only 1 populations, the tibble will not contain an "Overall" column.
#' P-values with "No contingency table" are replaced with "1" and p-values = "0" are replaced with 1 / (n_Batches x n_Iterations per batch).
#' 
#' @seealso 
#' [genepop::genepop-package()]
#' [genepop::test_LD()]
#' [GCLr::test_LD()]
#' 
#' @examples
#' \dontrun{
#' genepop_ld <- GCLr::read_genepop_dis(file = "~/R/test.txt.DIS", loci = loci, sillyvec = sillyvec)
#' }
#'  
#' @export
read_genepop_dis <- function(file, loci = NULL, sillyvec = NULL) {

  dis <- scan(file, what = '', sep = "\n")
  
  npops <- as.numeric(strsplit(dis[grep("Number of populations detected :", dis)], split = "Number of populations detected :")[[1]][2])
  
  nloci <- as.numeric(strsplit(dis[grep("Number of loci detected        :", dis)], split = "Number of loci detected        :")[[1]][2])

  batches <- as.numeric(strsplit(dis[grep("	Batches              :", dis)], split = "\tBatches              : ")[[1]][2])

  iterations <- as.numeric(strsplit(dis[grep("	Iterations per batch : ", dis)], split = "\tIterations per batch : ")[[1]][2])

  repzero <- format(1/(batches * iterations), scientific = FALSE, digits = 6)

  popstart <- grep("Pop             Locus#1  Locus#2    P-Value      S.E.     Switches", dis) + 2
  
  ncomps <- choose(nloci, 2)
  
  popend <- popstart + ncomps * npops - 1
  
  #Loci
  if (!is.null(loci)) {
    
    loc1 <- sapply(seq(length(loci)-1), function(i){
      
      loci[1:i]
    
    }) %>% unlist()
    
    loc2 <- sapply(seq(length(loci)), function(i){
      
      rep(loci[i], i-1)
    
    }) %>% unlist()
    
    pop_df0 <- tidyr::separate(data = tibble::tibble(dat = dis[popstart:popend]), col = dat, sep = "[[:blank:]]+", into = c("Pop", "Locus1", "Locus2", "PValue", NA, NA), remove = TRUE) %>% 
      dplyr::mutate(`Locus1`= rep(loc1, npops), `Locus2`= rep(loc2, npops), PValue = gsub(PValue, pattern = "No", replacement = "1")) %>% 
      dplyr::mutate(PValue = gsub(PValue, pattern = "$$0", replacement = repzero) %>% as.numeric())
    
  } else {
    
    pop_df0 <- tidyr::separate(data = tibble::tibble(dat = dis[popstart:popend]), col = dat, sep = "[[:blank:]]+", into = c("Pop", "Locus1", "Locus2", "PValue", NA, NA), remove = TRUE) %>% 
      dplyr::mutate(PValue = gsub(PValue, pattern = "No", replacement = "1")) %>% 
      dplyr::mutate(PValue = gsub(PValue, pattern = "$$0", replacement = repzero) %>% as.numeric())
    
  } # end Loci
  
  if (!is.null(sillyvec)) {
    
    pop_names <- lapply(sillyvec, function(silly){
      
      rep(silly, ncomps)
      
    }) %>% 
      unlist()
    
  } else {
    
    pop_names <- pop_df0 %>% 
      dplyr::pull(Pop)
    
  } # end sillyvec

  
  pop_df <- pop_df0 %>%
    dplyr::mutate(Pop = pop_names) %>% 
    dplyr::group_by(dplyr::across(-PValue)) %>%
    dplyr::mutate(row_id = 1:dplyr::n()) %>% 
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = Pop, values_from = PValue) %>% 
    dplyr::select(-row_id)
  
  # npops
  if (npops == 1) {
    
      summary_df <- pop_df
  
  } else {
    
    pop_df <- pop_df
    
    locstart <- popend + 6
    
    locend <- locstart + ncomps - 1
    
    loc_df <- tibble::tibble(dat = dis[locstart:locend]) %>%
      dplyr::mutate(dat = substr(dat, start = 46, stop = 54)) %>% 
      dplyr::mutate(dat_clean = dplyr::case_when(dat == "" ~ NA_character_,
                                                 substr(dat, 1, 1) == ">" ~ "0.000000",
                                                 TRUE ~ dat)) %>% 
      dplyr::mutate(dat_num = suppressWarnings(as.numeric(dat_clean)))

    summary_df = pop_df %>% 
      dplyr::mutate(Overall = loc_df %>% dplyr::pull(dat_num))
    
  } # end npops
    
     return(summary_df)

}