
# Find most likely sequence of IBD states in a single chromosome
.viterbiPath = function(.data, k1, a, Xchrom = FALSE, sex = NULL) {
  if(length(unique.default(.data$chrom)) > 1)
    stop2("`.viterbiPath()` expects a single chromosome")

  if(!all(c("g1", "g2", "cm", "freq1") %in% names(.data)))
    stop2("`.viterbiPath()` expects columns `g1`, `g2`, `cm` and `freq1`")

  n = nrow(.data)
  d = c(0, diff(.data$cm)) / 100

  emiss = as.matrix(.data[c("emission0","emission1")])

  # Viterbi algorithm -------------------------------------------------------

  inits = c(nonIBD = 1 - k1, IBD = k1)

  vit = ptr = matrix(NA_real_, 2, n)
  vit[,1] = log(inits) + log(emiss[1,])

  for(i in 2:n) {
    tm = log(trans(d[i], k1, a))
    for(s in 1:2) {
      scores = vit[, i-1] + tm[, s]
      ptr[s,i] = which.max(scores)
      vit[s,i] = max(scores) + log(emiss[i,s])
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
#' Identify genomic segments shared identical-by-descent (IBD) between two individuals
#' from SNP marker data. The underlying method uses a hidden Markov model (HMM) and
#' the Viterbi algorithm to infer the most likely IBD state (0 = non-IBD, 1 = IBD)
#' at each marker.
#'
#' @param data Data frame with required columns `chrom`, `cm`, `a1` and `freq1`.
#' @param ids Genotype columns (default: last 2 columns).
#' @param k1,a HMM parameters. See [fitHMM()] for how to estimate these.
#' @param prepped A logical indicating if the input data has undergone internal
#'   prepping. Can be ignored by most users.
#'
#' @returns Data frame with IBD segments.
#' @seealso [plotIBD()]
#'
#' @examples
#' findSegments(cousinsDemo, k1 = 0.2, a = 5)
#'
#' @export
findSegments = function(data, ids = NULL, k1, a, prepped = FALSE) {

  .data = if(prepped) data else prepForHMM(data, ids = ids)

  # X and sex
  Xchrom = attr(.data, "Xchrom")
  sex = attr(.data, "sex")

  seglist = list()
  chroms = unique(.data$chrom)

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
