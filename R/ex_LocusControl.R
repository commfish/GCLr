#' Example LocusControl
#' 
#' This is an example locus control object. 
#' 
#' @format A tibble data frame with 98 rows and 6 columns:
#' \describe{
#'   \item{MarkerSuite}{markersuite name from Loki}
#'   \item{locusnames}{Gene Conservation Lab locus names}
#'   \item{Publishedlocusnames}{published locus names}
#'   \item{nalleles}{number of alleles for each locus}
#'   \item{ploidy}{ploidy of each locus: 1 = haploid, 2 = diploid}
#'   \item{alleles}{list with an element for each locus. Each element is a 2 column tibble with allele (allele order) and call (allowable genotypes)}
#' }
#'
"ex_LocusControl"
