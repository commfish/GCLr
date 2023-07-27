#' Generate pooling tests
#' 
#' This function generates all possible pooling tests for a list of silly codes. The output is plugged into fishers_test(freq = x, loci = x, tests = combos) 
#' 
#' @param sillys  vector of "silly" elements.
#' 
#' @param min_length The minimum length of combinations to generate. Default is 2 for pairs of sillys.
#' 
#' @param max_length The maximum length of combinations to generate. 
#' 
#' @return a list item containing each of the tests
#' 
#' @examples
#' silly_elements <- c("silly1", "silly2", "silly3", "silly4", "silly5", "silly6")
#' min_length <- 2
#' max_length <- 6
#' combos <- fisher_combinations(silly_elements, min_length, max_length)
#' 
#' @export
fisher_combinations <- function(sillys, min_length = 2, max_length) {
  result <- sapply(min_length:max_length, function(n) {
    combs <- t(combn(sillys, n))
    apply(combs, 1, list) %>% unlist(recursive = FALSE)
  })
  unlist(result, recursive = FALSE)
}
