#' Dataset with autosomal SNP genotypes for two cousins
#'
#' Simulated genotypes for two individuals at the autosomal kinship SNPs from
#' the FORCE panel (Tillmar et al., 2021). The data was generated with the
#' `ibdsim2` package, assuming a relationship of first cousins once removed.
#'
#' @format A tibble with 3,915 rows and 9 variables:
#' * `CHROM`: Chromosome label
#' * `MARKER`: SNP identifier
#' * `MB`: Physical position in megabases
#' * `CM`: Map position in centiMorgans
#' * `A1`: First SNP allele
#' * `A2`: Second SNP allele
#' * `FREQ1`: Population frequency of `A1`
#' * `ID1`: Genotype of individual 1
#' * `ID2`: Genotype of individual 2
#'
#' @docType data
#' @keywords datasets
#'
#' @references Tillmar et al. *The FORCE Panel: An All-in-One SNP Marker Set for
#'   Confirming Investigative Genetic Genealogy Leads and for General Forensic
#'   Applications*. Genes. (2021)
#'
#' @examples
#' cousinsDemo
#'
"cousinsDemo"


#' Dataset with X-chromosomal SNP genotypes for two brothers
#'
#' Simulated genotypes for two brothers at the X-chromosomal SNPs included in
#' the FORCE panel (Tillmar et al., 2021). The data was generated with the
#' `ibdsim2` package.
#'
#' @format A tibble with 246 rows and 9 variables:
#' * `CHROM`: Chromosome label
#' * `MARKER`: SNP identifier
#' * `MB`: Physical position in megabases
#' * `CM`: Map position in centiMorgans
#' * `A1`: First SNP allele
#' * `A2`: Second SNP allele
#' * `FREQ1`: Population frequency of `A1`
#' * `ID1`: Genotype of individual 1
#' * `ID2`: Genotype of individual 2
#'
#' @docType data
#' @keywords datasets
#'
#' @references Tillmar et al. *The FORCE Panel: An All-in-One SNP Marker Set for
#'   Confirming Investigative Genetic Genealogy Leads and for General Forensic
#'   Applications*. Genes. (2021)

#' @examples
#' brothersX
#'
"brothersX"
