#' All-in-one workflow for finding IBD segments
#'
#' This function conveniently wraps the key steps of the package. It first fits
#' a continuous-time HMM to the data ([fitHMM()]), then identifies IBD segments
#' ([findSegments()]), and finally computes the marker-wise posterior IBD
#' probability at each marker locus ([ibdPosteriors()]). The result can be
#' passed straight to [plotIBD()] for visualisation.
#'
#' @param data Data frame with required columns `chrom`, `cm`, `a1` and `freq1`
#'   (case insensitive). Alternatively, a `ped` object, in which case the SNP
#'   data is extracted internally.
#' @param ids Character vector indicating genotype columns of `data` (default:
#'   last 2 columns).
#' @param k1,a HMM parameters passed on to [fitHMM()]. Supplying a value fixes
#'   the parameter; if NULL (default), the parameter is estimated.
#' @param err Error rate; a single number in `[0,1]` (default: 0).
#' @param thompson A logical passed on to [fitHMM()]. Default: FALSE.
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
findIBD = function(data, ids = NULL, k1 = NULL, a = NULL, thompson = FALSE,
                   verbose = TRUE)  {
  st = Sys.time()

  if(is.ped(data) || is.pedList(data)) {
    if(verbose)
      cat("Extracting SNP data from pedigree\n")
    data = getSNPdata(data)
  }

  .data = prepForHMM(data, ids = ids, err = err)

  if(verbose)
    cat("Individuals:", toString(attr(.data, "ids")), "\n")

  params = fitHMM(.data, k1 = k1, a = a, prepped = TRUE, thompson = thompson,
                  verbose = verbose)

  segs = findSegments(.data, k1 = params$k1, a = params$a, prepped = TRUE,
                      verbose = verbose)

  post = ibdPosteriors(.data, k1 = params$k1, a = params$a, prepped = TRUE,
                       verbose = verbose)

  if(verbose)
    cat("Analysis complete in", format(Sys.time() - st, digits = 3), "\n")

  list(ids = attr(.data, "ids"), k1 = params$k1, a = params$a, segments = segs,
       posteriors = post)
}
