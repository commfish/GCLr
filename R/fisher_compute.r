fisher_compute <- function(freq, loci, prec = 4){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function is a wrapper for stats::fisher.test() and runs a Fisher's Exact Test for homogeneity of allele frequencies.
  #   Null hypothesis: there are no significant differences in allele freqencies.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   freq - a tibble of allele freqencies produced by calc_freq_pop() for the collections you want to test.
  #
  #   loci - a character vector of locus names
  #   
  #   prec - the precision of the output pvalues (i.e., the number of significant digits)
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #    Produces a tibble of overall pvalues produced by combining individual locus pvalues using Fisher's method and a nested tibble of pvalues by locus.
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   load("V:/Analysis/2_Central/Chinook/Cook Inlet/2019/2019_UCI_Chinook_baseline_hap_data/2019_UCI_Chinook_baseline_hap_data.RData")
  #   old2new_locuscontrol()
  #   old2new_gcl(sillyvec = c("KKILL05","KKILL06"), save_old = FALSE)
  #
  #   freq <- calc_freq_pop(sillyvec = c("KKILL05","KKILL06"), loci = loci443)
  # 
  #   fisher_compute(freq = freq, loci = loci443, prec = 4)
  # 
  # Note~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function is called on by fishers_test to do multiple tests at once.
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	

  #Getting pvalues for each locus.
  pval <- sapply(loci, function(locus){
    
   my.freq <- freq %>% 
      filter(locus == !!locus) %>% 
     select(silly, allele_no, freq) 
   
   freq.mat <- matrix(my.freq$freq, ncol = max(my.freq$allele_no), dimnames = list(unique(my.freq$silly), unique(my.freq$allele_no)), byrow = TRUE)
   
   suppressWarnings(fisher.test(freq.mat, workspace = 8000000, hybrid = TRUE))$p.value
    
  })
  
  pval[pval==0] <- .Machine$double.xmin #Replacing all hard zeros with the smallest nonnegative number possible in R.
  
  #Using Fisher's method to get overall pvalue
  overall <- pchisq(q =- 2*sum(log(pval)), df = 2*length(loci), lower.tail = FALSE) %>% 
    round(prec)
  
  bylocus <- tibble::tibble(locus = names(pval), pval = pval %>% round(prec))
  
  test <- paste0(unique(freq$silly), collapse = ".")
  
  output <- tibble::tibble(test_sillys = test,  overall = !!overall, bylocus = list(!!bylocus))

  return(output)
  
}
 


