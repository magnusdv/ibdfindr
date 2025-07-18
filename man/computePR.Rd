% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computePR.R
\name{computePR}
\alias{computePR}
\title{Precision and Recall for IBD segment calls}
\usage{
computePR(call, truth, details = FALSE)
}
\arguments{
\item{call, truth}{Data frames with IBD segments, each with columns \code{chrom},
\code{startCM} and \code{endCM}.}

\item{details}{A logical indicating if additional details should be included
in the output.}
}
\value{
A data frame with columns \code{Precision} and \code{Recall}. If \code{details = TRUE}, additional columns \code{F1}, \code{TP} (true positives), \code{FP} (false positives)
and \code{FN} (false negatives) are included.
}
\description{
Computes the precision and recall of IBD segment calls (typically from
\code{\link[=findIBD]{findIBD()}}) against a truth set of IBD segments.
}
\examples{

# Built-in X example
ibd = findIBD(brothersX)

# True segments (see code in `data-raw/brothersX.R`)
truth = data.frame(chrom = 23,
                   startCM = c(0, 66.841, 138.834),
                   endCM = c(10.867, 120.835, 164.398))

computePR(ibd$segments, truth)
plotIBD(ibd, refSegs = truth)

}
