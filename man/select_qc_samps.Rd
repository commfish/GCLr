% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_qc_samps.R
\name{select_qc_samps}
\alias{select_qc_samps}
\title{Select High Quality QC Samples}
\usage{
select_qc_samps(project, username, password, ncores = 4)
}
\arguments{
\item{project}{the project name as it is spelled in Loki.}

\item{username}{Your user name for accessing LOKI.}

\item{password}{Your password for accessing LOKI.}

\item{ncores}{A numeric value for the number of cores to use in a \pkg{foreach} \verb{\%dopar\%} loop (default = 4).}
}
\value{
a tibble of all project samples containing the following variables: \verb{Container Type}, \code{WBID}, \code{Barcode}, \verb{Well Position}, \code{SILLY}, \verb{Fish ID}, \code{Tissue}, \code{SILLY_Fish_ID}, \code{select}
The function also writes out a csv file of the selection tibble (e.g., K215_QC_selection.csv) and another file of samples that were not removed for missing genotypes (e.g., K215_good_samples.csv) in case manual adjustments need to be made to the selection.
}
\description{
This function pulls project genotypes from Loki and removes samples with genotypes at less than 80\% of loci, then selects QC samples by plate (see details).
}
\details{
When selecting samples for QC analysis, the function tries to select the normal QC card positions first. (i.e., \code{c("A3", "B4", "C5", "D6", "E7", "F8", "G9", "H10")})
If any of those samples were removed for missing genotypes, the function selects a random sample from the same row as a replacement.  No samples are selected from extraction plates with less than 3 columns.
}
\examples{
\dontrun{

select_qc_samps(project = "K215", username = "awbarclay", password = "mypassword", ncores = parallel::detectCores())
 
}

}
