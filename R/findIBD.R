#' All-in-one workflow for finding IBD segments
#'
#' This function conveniently wraps the key steps of the package. It first fits
#' a continuous-time HMM to the data ([fitHMM()]), then identifies IBD segments
#' ([findSegments()]), and finally computes marker-wise posterior IBD
#' probabilities ([ibdPosteriors()]). The result can be passed straight to
#' plotIBD() for visualisation.
#'
#' @param data Data frame with required columns `chrom`, `cm`, `a1` and `freq1`
#'   (case insensitive).
#' @param ids Genotype columns (default: last 2 columns).
#' @param k1,a HMM parameters passed on to [fitHMM()]. Supplying a value fixes
#'   the parameter; if NULL (default), the parameter is estimated.
#' @param thompson A logical passed on to [fitHMM()].
#' @param verbose A logical, by default TRUE.
#'
#' @returns A list with the following elements:
#'   * `k1`: HMM parameter (estimated or provided)
#'   * `a`: HMM parameter (estimated or provided)
#'   * `segments`: Data frame with IBD segments.
#'   * `posteriors`: Data frame with posterior IBD probs at each marker.
#'
#' @seealso [fitHMM()], [findSegments()], [ibdPosteriors()], [plotIBD()]
#'
#' @examples
#' ibd = findIBD(brothersX)
#' plotIBD(ibd)
#'
#' @export
findIBD = function(data, ids = NULL, k1 = NULL, a = NULL, thompson = TRUE,
                   verbose = TRUE)  {

  .data = prepForHMM(data, ids = ids)

  params = fitHMM(.data, k1 = k1, a = a, prepped = TRUE, thompson = thompson,
                  verbose = verbose)

  if(verbose) {
    cat(sprintf("Fitting HMM done:\n  k1 = %.3f, a = %.3f\n", params$k1, params$a))
    cat("Finding IBD segments...")
  }

  segs = findSegments(.data, k1 = params$k1, a = params$a, prepped = TRUE)
  if(verbose) {
    cat("done\n ", nrow(segs), "segments found\n")
    cat("Calculating IBD posteriors...")
  }

  post = ibdPosteriors(.data, k1 = params$k1, a = params$a, prepped = TRUE)
  if(verbose) {
    cat("done\n")
  }

  list(ids = attr(.data, "ids"), k1 = params$k1, a = params$a,
       segments = segs, posteriors = post)
}
