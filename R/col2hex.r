#' Convert R color to hexadecimal
#' 
#' This function takes an R color name and convert it to hexidecimal format.
#' 
#' @param color an R color name
#' 
#' @import grDevices
#' @import dplyr
#' @import magrittr
#' @import tibble
#' @import tidyr
#' 
#' @return a tibble containing the following variables: 
#'          \itemize{
#'              \item \code{col}: R color name
#'              \item \code{hex}: hexadecimal color
#'              \item \code{r}: red saturation
#'              \item \code{g}: green saturation
#'              \item \code{b}: blue saturation
#'              }
#' 
#' @examples
#' col2hex(color = "red")
#' 
#' @export
#' 
col2hex <- function(color){
  
  c <- grDevices::col2rgb(color)
  
  cbind(col, sprintf("#%02X%02X%02X", c[1], c[2], c[3]), sprintf("%03d %03d %03d", c[1], c[2], c[3])) %>% 
    tibble::as_tibble() %>% 
    dplyr::rename(hex = V2) %>%
    tidyr::separate(V3, into = c("r", "g", "b"), sep = " ")
  
}
