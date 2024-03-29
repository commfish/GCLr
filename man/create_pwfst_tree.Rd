% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_pwfst_tree.r
\name{create_pwfst_tree}
\alias{create_pwfst_tree}
\title{Create Pairwise Fst Tree}
\usage{
create_pwfst_tree(
  sillyvec,
  loci,
  dir = NULL,
  nboots = 1000,
  ncores = 4,
  returnbootstrapFst = FALSE,
  LocusCtl = LocusControl
)
}
\arguments{
\item{sillyvec}{A character vector of silly codes without the ".gcl" extension.}

\item{loci}{A character vector of locus names.}

\item{dir}{Directory to save the \code{PairwiseFstTree} list object using \code{\link[base:dput]{base::dput()}}. If NULL (default), no file will be saved.}

\item{nboots}{A numeric value indicating the number of bootstrap iterations (default = 1000).}

\item{ncores}{A numeric value for the number of cores to use in a \pkg{foreach} \verb{\%dopar\%} loop (default = 4).
The number of cores cannot exceeds the number on your device (\code{\link[parallel:detectCores]{parallel::detectCores()}}).}

\item{returnbootstrapFst}{A logical value indicating whether to return the Fst matrix for each bootstrap iteration (default = \code{FALSE}).}

\item{LocusCtl}{an object created by \code{\link[=create_locuscontrol]{create_locuscontrol()}} (default = LocusControl).}
}
\value{
A list with the following 4 or 5 components:
\itemize{
\item \code{tree}: a neighbor joining tree list of 4 created by \code{\link[ape:nj]{ape::nj()}} containing:
\itemize{
\item \code{edge}: numeric matrix
\item \code{edge.length}: numeric vector
\item \code{tip.lable}: character vector of population names, inherited from \code{sillyvec}
\item \code{Nnode}: integer
}
\item \code{bootstrap}: numeric vector of node boodstrap values
\item \code{PairwiseFst}: numeric matrix (\code{length(sillyvec)} x \code{length(sillyvec)}) of pairwise Fst values
\item \code{vc}: list of \code{choose(n = length(sillyvec), k = 2)} with pairwise variance components from \code{\link[hierfstat:varcomp]{hierfstat::varcomp()}}
\item \code{BootstrapFst}: list of \code{nboots} with numeric matrix of pairwise Fst values for each bootstrap iteration
(optional, only returned if \code{returnbootstrapFst = TRUE})
}
This list object is saved in \code{dir} using \code{\link[base:dput]{base::dput()}}, and is named \code{paste(dir,"\\\\", length(sillyvec), "Pops", length(loci), "Loci_", "PairwiseFstTree.txt", sep = "")}
}
\description{
This function generates a matrix of pairwise Fst values, a neighbor joining tree, bootstrap values for tree nodes, and variance components.
It can utilize multicore processing with \pkg{foreach} to speed up the calculation of variance components.
}
\details{
Older versions of this function used to write out an FSTAT .dat file, but this is no longer the case.
The function now calls on \code{\link[=create_hierfstat_data]{create_hierfstat_data()}} to create a \pkg{hierfstat} data object internally.
Depending on the size of your baseline and number of bootstrap iterations, this function can take a while to run.
}
\examples{
\dontrun{
  
  sillyvec <- GCLr::base2gcl(GCLr::ex_baseline)
  
  groupvec <- GCLr::ex_baseline \%>\%
    dplyr::group_by(collection) \%>\%
    dplyr::summarize(group = unique(repunit)) \%>\%
    dplyr::mutate(group = factor(group) \%>\% as.numeric()) \%>\%
    dplyr::pull(group)
  
  loci <- GCLr::ex_baseline[,-c(1:5)] \%>\%
    names() \%>\%
    gsub(pattern = "*\\\\.1", x = ., replacement = "") \%>\%
    unique()
  
  GCLr::create_pwfst_tree(sillyvec = sillyvec, loci = loci, dir = NULL, nboots = 1000, ncores = parallel::detectCores(), returnbootstrapFst = FALSE)
  
}

}
