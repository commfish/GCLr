#' Summarize LD Results
#'
#' This function takes the output from [GCLr::test_LD()], filters tests by the supplied significance level (alpha),
#' and makes a bar plot of the proportion of pops with significant test results for each locus pair.
#'
#' @param LDresults The output from [GCLr::test_LD()]; a tibble with 3 variables:
#'                  \itemize{
#'                     \item \code{Locus_pair}
#'                     \item \code{Pop}
#'                     \item \code{Pvalue}
#'                     }
#'                 
#' @param alpha A numeric value indicating the significance level for filtering the [GCLr::test_LD()] output; default is 0.05.
#' 
#' @param prop_sign_pops A numeric value indicating the proportion of populations with significant test results; default is 0.5.
#'                       Only locus pairs that have > prop_sign_pops will be plotted.
#'
#' @return An interactive bar plot with the proportion of populations with significant tests on the y-axis and the locus pairs on the x-axis.
#' 
#' @seealso [GCLr::gcl2genepop()]
#' @seealso [GCLr::test_LD()]
#'
#' @examples
#' \dontrun{
#' 
#' sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
#' 
#' loci <- GCLr::ex_LocusControl$locusnames[-c(10, 12, 13, 32, 33, 97, 98)]
#' 
#' dir.create(path = "~/GENEPOP")
#' 
#' GCLr::gcl2genepop(sillyvec = sillyvec, loci = loci, path = path.expand("~/GENEPOP"), VialNums = TRUE, usat = FALSE, 
#'                   ncores = parallel::detectCores(), npops = 1, LocusCtl = GCLr::ex_LocusControl) 
#' 
#' LDresults <- GCLr::test_LD(genepopFiles = list.files(path.expand("~/GENEPOP")), path = path.expand("~/GENEPOP"), batches = 14, iterations = 1, ncores = parallel::detectCores())
#'
#' GCLr::summarize_LD(LDresults, alpha = 0.05, prop_sign_pops = 0.5)
#' 
#' }
#'
#' @export
summarize_LD <- function(LDresults, alpha = 0.05, prop_sign_pops = 0.5){
 
  npops <- length(LDresults$Pop %>% unique())

  PLD <- LDresults %>% 
    dplyr::filter(Pvalue < alpha) %>% 
    dplyr::group_by(Locus_pair) %>% 
    dplyr::summarize(prop_pops_LD = length(Pop)/npops) %>% 
    dplyr::filter(prop_pops_LD > prop_sign_pops)
  
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
    dplyr::rename(`Locus Pair` = Locus_pair, `Proportion of Populations` = prop_pops_LD)
  
  plotly::ggplotly(plot_df %>% 
                     ggplot2::ggplot(ggplot2::aes(x = `Locus Pair`, y = `Proportion of Populations`))+
                     ggplot2::geom_bar(stat = "identity")+
                     ggplot2::ylim(0, 1)+
                     ggplot2::ggtitle(label = my.text))
  
}