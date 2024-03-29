% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_extraction_info.R
\name{get_extraction_info}
\alias{get_extraction_info}
\title{Get DNA Extraction Information}
\usage{
get_extraction_info(
  plate_ids = NULL,
  sillyvec = NULL,
  username,
  password,
  file = NULL
)
}
\arguments{
\item{plate_ids}{a numeric vector of plate IDs (e.g. plate_ids <- c(4142, 7221, 7222, 8589)) (default = NULL)}

\item{sillyvec}{optional, a character vector of silly codes (e.g. sillyvec <- c("KQUART06","KQUART08","KQUART09")), used to get plate_ids if you do not know them (default = NULL)}

\item{password}{\itemize{
\item your password used to access LOKI - contact Eric Lardizabal if you don't have a password for LOKI
}}

\item{file}{the file path, including .csv extension, for writing out a csv file of the output. (default = NULL)}
}
\value{
The function outputs a tibble containing the following 5 variables:
\itemize{
\item \code{WELL_NO}:
\item \code{FK_PLATE_ID}: the extraction plate ID
\item \code{SILLY_CODE}: the collection silly code
\item \code{FISH_NO}: the collection fish number
\item \code{TISSUETYPE}: e.g., Axillary Process, Fin, Heart, etc.
}
}
\description{
This function connects to LOKI and pulls the well information for each DNA extraction plate.
}
\details{
This function can pull extraction information from Loki using either the plate IDs (if known) or using silly codes (if plate IDs aren't known).
This function was designed to pull all wells for each plate. If you supply \code{sillyvec}, it will first query LOKI to figure out which plate_ids have those sillys, and then it will query LOKI again for all wells in those plate_ids.
When using plate IDs, the output will only contain extraction information for the supplied plate IDs.
If using silly codes, the output will contain information for all extraction plates for the supplied silly codes.
}
\examples{
\dontrun{
 
password = "************"

#By plate ID
GCLr::get_extraction_info(plate_ids = c(4142, 7221, 7222, 8589), sillyvec = NULL, username = "krshedd", password = password, file = path.expand("~/plateID_extraction_info.csv") )

#By silly code
GCLr::get_extraction_info(plate_ids = NULL, sillyvec = c("KKENM06", "KQUART06", "KCRESC06", "KRUSSR08", "KQUART08", "KDAVE08", "KSLIK08", "KKENU09", "KQUART09", "KKENN09"), username = "krshedd", password = password, file = path.expand("~/silly_extraction_info.csv") )

}
}
