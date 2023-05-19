read_genepop_dis <- function(file, loci = NULL, sillyvec = NULL) {
#################################################################################################################
#
# This function reads in the output from a GENEPOP disequilibrium ("*.dis") file and returns a tibble with p-value 
# columns for each population and a column of overall population p-values. If the disequilibrium file contains tests 
# for only 1 population, the tibble will not contain on 'Overall' column.
#
# P-values with "No contingency table" are replaced with "1" and p-values = "0" are replaced with 1/(#Batches*#Iterations per batch).
# 
# Input parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# file <- "V:/Analysis/2_Central/Coho/Cook Inlet/2018/Baseline/Genepop/CI94pops82loci.DIS"
#   ~ the full file path, including the ".dis" extension.  Make sure the file has not been modified.
#
# attach("V:/Analysis/2_Central/Coho/Cook Inlet/2018/Baseline/CI_Coho2018Baseline.RData")
# 
# loci <- loci82
#   ~ default is NULL or a vector of locus names the same order as the genepop input file used for the tests. 
# If a vector of locus names is supplied the Locus1 and Locus2 columns will contain the correct locus names
# otherwise the locus names will be the same as in the .DIS file (truncated to 8 characters)
#
# sillyvec <- sillyvec94
#   ~ default is NULL or a vector of population sillys the same order as in the genepop input file used for the tests. 
#     If a vector of silly's is supplied the header names for each pop in the output tibble will be the correct sillys, 
#     otherwise, if no vector is supplied, the header names will be the SillySource of the last fish in the population.
#
# detach("file:V:/Analysis/2_Central/Coho/Cook Inlet/2018/Baseline/CI_Coho2018Baseline.RData")
# Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# tibble
# nrow = number of pairwise combinations of loci (i.e. choose(n = nloci, k = 2)).
# ncol = number of populations + 3
#   ~ The first two columns are the locus names for the pair (Locus1 and Locus2)
#   ~ Subsequent columns contain the LD p-values for that locus pair for a given population
#   ~ The last column contains the overall p-values for all populations. 
#     If only one population is tested, the tibble will not include an 'Overall' column. 
#
# Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How many populations have a p-value of < 0.05 for a given locus pair?
# Create an additional column with the number of pops below a given p-value threshold (i.e. 0.05)
# LD <- read_genepop_dis(file = "V:/Analysis/3_AYK/Chum/WASC/Genepop/WASCchum.DIS",loci = NULL)
# LD$npopsfail <- apply(LD[ , 3:dim(LD)[2]] < 0.05, 1, function(locuspair) {sum(locuspair)} )
#
#################################################################################################################

  require(tidyverse)

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
  if(!is.null(loci)) {
    
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
    
  } #end Loci
  
  if(!is.null(sillyvec)){
    
    pop_names <- lapply(sillyvec, function(silly){
      
      rep(silly, ncomps)
      
    }) %>% 
      unlist()
    
  }else{
    
    pop_names <- pop_df0 %>% 
      dplyr::pull(Pop)
    
  } #end sillyvec

  
  pop_df <- pop_df0 %>%
    dplyr::mutate(Pop = pop_names) %>% 
    dplyr::group_by(dplyr::across(-PValue)) %>%
    dplyr::mutate(row_id = 1:n()) %>% 
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = Pop, values_from = PValue) %>% 
    dplyr::select(-row_id)
  
  #npops
  if(npops==1) {
    
      summary_df <- pop_df
  
  } else {
    
    pop_df <- pop_df
    
    locstart <- popend + 6
    
    locend <- locstart + ncomps - 1
    
    loc_df <- tibble::tibble(dat = dis[locstart:locend]) %>%
      dplyr::mutate(dat = substr(dat, start = 46, stop = 54)) %>% 
      dplyr::mutate(dat = gsub(pattern = "Highly si", x = dat, replacement = "0.000000") %>% as.numeric()) 
    
    summary_df = pop_df %>% 
      dplyr::mutate(Overall = loc_df %>% dplyr::pull(dat))
    
  } #end npops
    
     return(summary_df)

}