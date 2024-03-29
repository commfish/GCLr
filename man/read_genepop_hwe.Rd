% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_genepop_hwe.r
\name{read_genepop_hwe}
\alias{read_genepop_hwe}
\title{Read Genepop HWE Output}
\usage{
read_genepop_hwe(file, sillyvec = NULL, summaryAsNumeric = TRUE)
}
\arguments{
\item{file}{The full file path to the \code{genepop} HWE output, including the ".P" extension.}

\item{sillyvec}{An optional character vector of silly codes without the ".gcl" extension (default = \code{NULL}).
If \code{NULL}, \code{Pop} names will come directly from the ("*.P") file and likely include "_fishID" extensions.
If supplying \code{sillyvec}, make sure it is the same \code{sillyvec} used in \code{\link[=gcl2genepop]{gcl2genepop()}}.}

\item{summaryAsNumeric}{A logical vector of length 1 indicating whether character p-values such as "High. sign." or " < "
should be coerced into numeric (default = \code{TRUE}).}
}
\value{
A list with 2 components:
\itemize{
\item \code{DataByPop}: a tibble with 7 columns containing the full HWE output:
\itemize{
\item \code{Pop}: silly code
\item \code{Locus}: locus name
\item \code{PValue}: HWE p-value
\item \code{SE}: standard error of the estimated p-value
\item \code{WC Fis}: Weir and Cockerham Fis estimate
\item \code{RH Fis}: Robertson and Hill Fis estimate
\item \code{Steps}: number of genotypic matrices considered (exact) or number of switches (MCMC)
}
\item \code{SummaryPValues}: a tibble containing p-values for each locus (row) and pop (column) including the
overall p-value across pops and loci (Fisher's method)
}
}
\description{
This function reads in output from a \code{genepop} Hardy-Weinberg test ("*.P") file.
}
\details{
In \code{genepop} version 4.6 or below, \code{summaryAsNumeric = TRUE}, forces character output such as "High. sign." and to "0.000".
In \code{genepop} version 4.7 or avove, \code{summaryAsNumeric = TRUE}, forces character output such as " < 0.1555" and to "0.1555".
If MCMC, replace hard 0 p-values with \code{1 / (batches * iterations)} and calculate overall pops and overall loci via Fisher's method (ChiSqaure).
if Exact test, p-values overall pops or overall loci are pulled directly from \code{genepop} ("*.P") file
}
\examples{
sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)

loci <- GCLr::ex_LocusControl$locusnames[-c(10, 12, 13, 32, 33, 97, 98)]

dir.create(path = "~/GENEPOP")

GCLr::gcl2genepop(sillyvec = sillyvec, loci = loci, path = path.expand("~/GENEPOP/my_genepop.txt"), VialNums = TRUE, usat = FALSE,
                  ncores = parallel::detectCores(), npops = NULL, LocusCtl = GCLr::ex_LocusControl)

genepop::test_HW(inputFile = path.expand("~/GENEPOP/my_genepop.txt"))

GCLr::read_genepop_hwe(file =  path.expand("~/GENEPOP/my_genepop.txt.P"), sillyvec = sillyvec)

}
\seealso{
\code{\link[genepop:genepop-package]{genepop::genepop-package()}}
\code{\link[genepop:Hardy-Weinberg]{genepop::test_HW()}}
\code{\link[=gcl2genepop]{gcl2genepop()}}
}
