#' Convert R color to hexadecimal
#' 
#' This function takes an R color name and convert it to hexidecimal format.
#' 
#' @param color an R color name
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
col2hex <- function(color){
  
  myc <- grDevices::col2rgb(color)
  
  cbind(col = color, sprintf("#%02X%02X%02X", myc[1], myc[2], myc[3]), sprintf("%03d %03d %03d", myc[1], myc[2], myc[3])) %>% 
    tibble::as_tibble(.name_repair = "universal") %>% 
    dplyr::rename(hex = ...2) %>%
    tidyr::separate(...3, into = c("r", "g", "b"), sep = " ")
  
}