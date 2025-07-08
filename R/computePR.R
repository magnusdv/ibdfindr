#' Precision and Recall for IBD segment calls
#'
#' Computes the precision and recall of IBD segment calls (typically from
#' [findIBD()]) against a truth set of IBD segments.
#'
#' @param call,truth Data frames with IBD segments, each with columns `chrom`,
#'   `startCM` and `endCM`.
#' @param details A logical indicating if additional details should be included
#'   in the output.
#'
#' @returns A data frame with columns `Precision` and `Recall`. If `details =
#' TRUE`, additional columns `F1`, `TP` (true positives), `FP` (false positives)
#' and `FN` (false negatives) are included.
#'
#' @examples
#'
#' # Built-in X example
#' ibd = findIBD(brothersX)
#'
#' # True segments (see code in `data-raw/brothersX.R`)
#' truth = data.frame(chrom = 23,
#'                    startCM = c(0, 66.841, 138.834),
#'                    endCM = c(10.867, 120.835, 164.398))
#'
#' computePR(ibd$segments, truth)
#' plotIBD(ibd, refSegs = truth)
#'
#' @export
computePR = function(call, truth, details = FALSE) {
  call = as.data.frame(call)
  truth = as.data.frame(truth)

  chroms = intersect(truth$chrom, call$chrom)
  cols = c("startCM", "endCM") # segment start / end

  TP = 0

  for(chr in chroms) {
    t1 = truth[truth$chrom == chr, cols, drop = FALSE]
    c1 = call[call$chrom == chr, cols, drop = FALSE]
    ov = rangeIntersect(list(as.matrix(t1), as.matrix(c1)))
    len = if(nrow(ov) == 0) 0 else sum(ov[,2] - ov[,1])
    TP = TP + len
  }

  totalT = sum(truth$endCM - truth$startCM)
  totalC = sum(call$endCM - call$startCM)
  precis = if(totalC > 0) TP / totalC else 0
  recall = TP / totalT

  res = data.frame(Precision = precis, Recall = recall)

  if(details) {
    res$F1 = 2 * precis * recall / (precis + recall)
    res$CallTotal = totalC
    res$TruthTotal = totalT
  }
  res
}

# Fast sweep-algorithm for intersection of two sorted interval lists
rangeIntersect = function(x) {
  x = lapply(x, as.matrix)
  nrws = lengths(x, use.names=FALSE)/2
  if (any(nrws == 0)) return(matrix(numeric(0), ncol=2))

  n = length(nrws)
  vec = unlist(x, recursive = TRUE, use.names = FALSE)
  startstop = rep.int(rep.int(c(1,-1), n), rep.int(nrws, rep.int(2,n)))

  ord = order(vec, -startstop, na.last=TRUE, decreasing=FALSE)
  m = vec[ord]
  startstop.ord = startstop[ord]

  ind.start = seq_along(vec)[cumsum(startstop.ord) == n]
  ind.end = ind.start + 1
  cbind(m[ind.start], m[ind.end], deparse.level=0)
}
