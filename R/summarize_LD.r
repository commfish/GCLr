summarize_LD <- function(LDresults, alpha = 0.05, prop_sign_pops = 0.5){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  #   This function takes the output from test_LD, filters tests by the supplied significance level (alpha), 
  #   and makes a bar plot of the proportion of pops with significant test results for each locus pair.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   LDresults - the output from test_LD; a tibble with 3 variables: 
  #               1) Locus_pair, 
  #               2) Pop, 
  #               3) Pvalue
  #
  #   alpha - a numeric string indicating the significance level for filtering the Test_LD.gcl output; default is 0.05
  #   
  #   prop_sign_pops - a numeric string indicating the proportion of populations with significant test results; default is 0.5. 
  #                    Only locus pairs that have > prop_sign_pops will be plotted. 
  # 
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  # An interactive bar plot with the proportion of populations with significant tests on the y-axis and the locus pairs on the x-axis.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #  attach("V:/Analysis/5_Coastwide/Chum/Baseline for South Peninsula MSA/South_AK_Peninsula_chum_baseline analysis.Rdata")
  #
  #  LDresults <- LD
  #
  #  detach()
  #
  #  summarize_LD(LDresults, alpha = 0.05, prop_sign_pops = 0.5)
  #   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  npops <- length(LDresults$Pop %>% unique())

  PLD <- LDresults %>% 
    filter(Pvalue < alpha) %>% 
    group_by(Locus_pair) %>% 
    summarize(prop_pops_LD = length(Pop)/npops) %>% 
    filter(prop_pops_LD > prop_sign_pops)
  
  if(dim(PLD)[1] == 0){
    
    my.text <- paste0("No locus pairs had significant tests (alpha = ", alpha, ") for ", prop_sign_pops*100, "% of pops.")
    
  }
  
  if(dim(PLD)[1] > 1){
    
    my.text <- paste0(dim(PLD)[1], " locus pairs had significant tests (alpha = ", alpha, ") for ", prop_sign_pops*100, "% of pops.")
    
    
  }
    
  if(dim(PLD)[1] == 1){
    
    my.text <- paste0(dim(PLD)[1], " locus pair had significant tests (alpha = ", alpha, ") for ", prop_sign_pops*100, "% of pops.")
    
    
  }
  
  plot_df <- PLD %>% 
    rename(`Locus Pair` = Locus_pair, `Proportion of Populations` = prop_pops_LD)
  
  plotly::ggplotly(plot_df %>% 
                     ggplot(aes(x = `Locus Pair`, y = `Proportion of Populations`))+
                     geom_bar(stat = "identity")+
                     ylim(0, 1)+
                     ggtitle(label = my.text))
  
}




