% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_gtscore_loki_import.R
\name{split_gtscore_loki_import}
\alias{split_gtscore_loki_import}
\title{Split a GTscore Import File into Multiple Files}
\usage{
split_gtscore_loki_import(file, nlines = 3e+05)
}
\arguments{
\item{file}{the file path to the GTscore file including .csv extension}

\item{nlines}{the maximum number of lines you want in each file (default = 3e+05)}
}
\description{
This function takes a GTscore import file and splits it up into
multiple files so genotypes can be imported into Loki in smaller batches (see details)
}
\details{
The function reads in a GTscore file and splits it up into multiple files with a maximum of number of lines = \code{nlines}
and writes them out in the same folder as the original file with letters at the end of the file names. For example,
if you have a file named S246_LOKI_input_all.csv with 1 million lines and you want to break it up into files with 300,000 lines each, the function will produce
3 files with 300,000 lines and 1 file with 100,000 lines, and name them "S246_LOKI_input_all_a.csv", "S246_LOKI_input_all_b.csv", "S246_LOKI_input_all_c.csv", and "S246_LOKI_input_all_d.csv."
}
\examples{
\dontrun{

GCLr::split_gtscore_loki_import(file = "C:/Users/awbarclay/Desktop/GTscore import files split/S246a_LOKI_input_all.csv", nlines = 3e+05) 

}

}
