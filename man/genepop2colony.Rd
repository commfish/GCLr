% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genepop2colony.r
\name{genepop2colony}
\alias{genepop2colony}
\title{Convert GENEPOP File to COLONY Input File}
\usage{
genepop2colony(path_genepop, path_colony)
}
\arguments{
\item{path_genepop}{The full file path of the GENEPOP file with either .gen or .txt extension.}

\item{path_colony}{The full file path of the COLONY file including .txt extension.}
}
\value{
A .txt file of genotypes for COLONY. The first column is the ID, followed by a 2 column genotypes (2 columns per locus).
}
\description{
Create COLONY genotype input file from a GENEPOP file.
}
\details{
This function is designed to convert a GENEPOP file to a COLONY genotype file.
If you are doing parentage analysis, you'll need to run this separately for parents and offspring.
If you are using adult sex data, you'll run the parents separately for Mum and Dad.
}
\examples{
GCLr::genepop2colony(path_genepop = system.file("genepop", "ex_genepop.txt", package = "GCLr"), path_colony = path.expand("~/example_colony_file.txt"))
}
