% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ibdPosteriors.R
\name{ibdPosteriors}
\alias{ibdPosteriors}
\title{IBD posteriors}
\usage{
ibdPosteriors(
  data,
  ids = NULL,
  k1,
  a,
  err = 0,
  prepped = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{data}{Data frame with required columns \code{chrom}, \code{cm}, \code{a1} and \code{freq1}.}

\item{ids}{Genotype columns (default: last 2 columns).}

\item{k1, a}{HMM parameters. See \code{\link[=fitHMM]{fitHMM()}} for how to estimate these.}

\item{err}{Error rate; a single number in \verb{[0,1]} (default: 0).}

\item{prepped}{A logical indicating if the input data has been internally
processed. Can be ignored by most users.}

\item{verbose}{A logical.}
}
\value{
Data frame similar to \code{data}, with a column \code{post} containing the
posterior IBD probability at each marker locus.
}
\description{
Computes the posterior probability of identity‐by‐descent (IBD) at each
marker locus via the HMM forward–backward algorithm.
}
\examples{
ibdPosteriors(cousinsDemo, k1 = 0.2, a = 5)

}
\seealso{
\code{\link[=plotIBD]{plotIBD()}}
}
