% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_proj_geno.R
\name{format_proj_geno}
\alias{format_proj_geno}
\title{Format project genotypes report}
\usage{
format_proj_geno(report = NULL, in.file = NULL, out.file = NULL)
}
\arguments{
\item{report}{a project genotypes report object created by \code{\link[=get_geno]{get_geno()}} (see details)}

\item{in.file}{the file path to an unmodified genotypes report .csv file created by \code{\link[=get_geno]{get_geno()}} (see details)}

\item{out.file}{the file path for the .xlsx output file with file name and extension.
If not supplied, the file named "GenotypeReportQC.xlsx" will be written to your current working directory (see \code{\link[base:getwd]{base::getwd()}})}
}
\value{
This function writes out a formatted genotypes report Excel file with three worksheets
\itemize{
\item \code{ProjectName_Genotype_Table}: the report produced by \code{\link[=get_geno]{get_geno()}} when pulling by project_name with "Count of 0/0" and "Success Rate" columns added and conditional formatting.
\item \code{Calculations}: number of fish with genotypes at less than 80\% of loci, the total number of project fish, and the proportion of fish with genotypes at greater than 80\% of loci.
\item \code{Avg success by plate}: the average success rate by plate ID
}
}
\description{
This function is used to format the gentoypes report produced by \code{\link[=get_geno]{get_geno()}}
}
\details{
This function requires either a report object or a report .csv file created by \code{\link[=get_geno]{get_geno()}} when \code{project_name} is supplied.
If a file and a report object are both supplied, the function will stop and throw and error message.
If \code{project_name} was not used to produce the file/object, the function will stop and throw a error message.
}
\examples{
\dontrun{

Geno_proj <- get_geno(project_name = "K205", file = "TestGenotypesReport_proj.csv", username = username, password = password, open.file = TRUE)

format_proj_geno(report = Geno_proj, in.file = NULL, out.file = NULL) #supplying a report object

format_proj_geno(report = NULL, in.file = "TestGenotypesReport_proj.csv", out.file = NULL) #supplying a report file

}

}
\seealso{
\code{\link[=get_geno]{get_geno()}}
}
