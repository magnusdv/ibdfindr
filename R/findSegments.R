
# Find most likely sequence of IBD states in a single chromosome
.viterbiPath = function(.data, k1, a, Xchrom = FALSE, sex = NULL) {
  if(length(unique.default(.data$chrom)) > 1)
    stop2("`.viterbiPath()` expects a single chromosome")

  if(!all(c("g1", "g2", "cm", "freq1") %in% names(.data)))
    stop2("`.viterbiPath()` expects columns `g1`, `g2`, `cm` and `freq1`")

  n = nrow(.data)
  d = c(0, diff(.data$cm)) / 100

  # Emission probabilities as 2*n matrix
  emProb = emissionMat(.data$freq1, g1 = .data$g1, g2 = .data$g2, Xchrom = Xchrom, sex = sex)

  # Viterbi algorithm -------------------------------------------------------

  inits = c(nonIBD = 1 - k1, IBD = k1)

  vit = ptr = matrix(NA_real_, 2, n)
  vit[,1] = log(inits) + log(emProb[,1])

  for(i in 2:n) {
    tm = log(trans(d[i], k1, a))
    for(s in 1:2) {
      scores = vit[, i-1] + tm[, s]
      ptr[s,i] = which.max(scores)
      vit[s,i] = max(scores) + log(emProb[s, i])
    }
  }

  # Traceback: Most likely path
  path = integer(n)
  path[n] = which.max(vit[, n])
  for(i in (n-1):1)
    path[i] = ptr[path[i+1], i+1]

  # Convert path to IBD states 0/1
  path - 1
}


#' Identify IBD segments
#'
#' @param data Data frame with required columns `chrom`, `cm`, `a1` and `freq1`.
#' @param ids Genotype columns (default: last 2 columns).
#' @param k1,a HMM parameters.
#'
#' @returns Data frame with IBD segments.
#' @seealso [plotIBD()]
#'
#' @examples
#' findSegments(cousinsDemo, k1 = 0.2, a = 5)
#'
#' @export
findSegments = function(data, ids = NULL, k1, a) {
  .data = prepForHMM(data, ids = ids)

  chroms = unique(.data$chrom)
  Xchrom = length(chroms) == 1 && chroms == 23
  sex = if(Xchrom) getsex(.data) else NULL

  seglist = list()

  # Loop through chroms
  for(chr in chroms) {
    subdat = .data[.data$chrom == chr, , drop  = FALSE]
    vpath = .viterbiPath(subdat, k1, a, Xchrom = Xchrom, sex = sex)

    # Extract IBD segments if any
    runs = rle(vpath)
    ibd = which(runs$values == 1)
    if(!length(ibd))
      next

    stops = cumsum(runs$lengths)
    starts = c(1, utils::head(stops, -1) + 1)

    segs = data.frame(chrom = chr,
                       start = subdat$cm[starts[ibd]],
                       end  = subdat$cm[stops[ibd]],
                       n     = runs$lengths[ibd])
    seglist[[chr]] = segs
  }

  res = do.call(rbind, seglist)
  rownames(res) = NULL
  res
}
