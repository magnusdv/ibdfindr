% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cousinsDemo}
\alias{cousinsDemo}
\title{Dataset with autosomal SNP genotypes for two cousins}
\format{
A tibble with 3,915 rows and 9 variables:
\itemize{
\item \code{CHROM}: Chromosome label
\item \code{MARKER}: SNP identifier
\item \code{MB}: Physical position in megabases
\item \code{CM}: Map position in centiMorgan
\item \code{A1}: First SNP allele
\item \code{A2}: Second SNP allele
\item \code{FREQ1}: Population frequency of \code{A1}
\item \code{ID1}: Genotype of individual 1
\item \code{ID2}: Genotype of individual 2
}
}
\usage{
cousinsDemo
}
\description{
Simulated genotypes for two individuals at the autosomal kinship SNPs from
the FORCE panel (Tillmar et al., 2021). The data was generated with the
\code{ibdsim2} package, assuming a relationship of first cousins.
}
\examples{
cousinsDemo

\dontrun{
r = findIBD(cousinsDemo)
plotIBD(r)
}

}
\references{
Tillmar et al. \emph{The FORCE Panel: An All-in-One SNP Marker Set for
Confirming Investigative Genetic Genealogy Leads and for General Forensic
Applications}. Genes. (2021)
}
\keyword{datasets}
