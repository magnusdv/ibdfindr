#' Example dataset with autosomal SNPs for two related individuals
#'
#' A dataset with 3,915 autosomal SNPs (from the FORCE panel) typed in two
#' related individuals.
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
#' @examples
#' cousinsDemo
"cousinsDemo"
