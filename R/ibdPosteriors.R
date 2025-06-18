# Compute posterior IBD probabilities for a single chromosome
.ibdPosteriors = function(.data, k1, a, Xchrom = FALSE, sex = NULL) {
  if(length(unique.default(.data$chrom)) > 1)
    stop("Function `.ibdPosterior()` expects a single chromosome.")

  n = nrow(.data)
  d = c(0, diff(.data$cm)) / 100

  # Matrix (2*n) of emission probabilities
  fr1 = .data$freq1
  names(fr1) = .data$marker %||% paste0("M", seq_len(n))
  emProb = emissionMat(fr1, g1 = .data$g1, g2 = .data$g2, Xchrom = Xchrom, sex = sex)

  # Forward-backward algorithm ----------------------------------------------

  # Initial state probabilities
  inits = c(nonIBD = 1 - k1, IBD = k1)

  # Forward (with logliks)
  fwd = matrix(0, 2, n)
  fwd[, 1] = (inits * emProb[, 1]) |> nrm()

  for(i in 2:n) {
    tm = trans(d[i], k1, a)
    fwd[ ,i] = ((tm %*% fwd[,i-1]) * emProb[,i]) |> nrm()
  }

  # Backward
  bwd = matrix(1, 2, n)
  for(i in (n-1):1) {
    tm = trans(d[i+1], k1, a)
    bwd[,i] = (tm %*% (emProb[,i+1] * bwd[,i+1])) |> nrm()
  }

  # Posterior IBD probs
  post = fwd * bwd
  post[2,]/colSums(post)
}


#' IBD posteriors
#'
#' Computes the posterior probability of identity‐by‐descent (IBD) at each
#' marker locus via the HMM forward–backward algorithm.
#'
#' @param data Data frame with required columns `chrom`, `cm`, `a1` and `freq1`.
#' @param ids Genotype columns (default: last 2 columns).
#' @param k1,a HMM parameters. See [fitHMM()] for how to estimate these.
#'
#' @returns Data frame similar to `data`, with a column `post` containing the
#'   posterior IBD probability at each marker locus.
#'
#' @seealso [plotIBD()]
#'
#' @examples
#' ibdPosteriors(cousinsDemo, k1 = 0.2, a = 5)
#'
#' @export
ibdPosteriors = function(data, ids = NULL, k1, a) {
  .data = prepForHMM(data, ids = ids)

  chroms = unique.default(.data$chrom)
  Xchrom = length(chroms) == 1 && chroms == 23
  sex = if(Xchrom) getsex(.data) else NULL

  # Initialize output objects
  dflist = list()

  # Loop through chroms
  for(chr in chroms) {
    subdat = .data[.data$chrom == chr, , drop  = FALSE]
    post = .ibdPosteriors(subdat, k1, a, Xchrom = Xchrom, sex = sex)

    # Add column with posteriors
    subdat$post = post
    dflist[[chr]] = subdat
  }

  res = do.call(rbind, dflist)
  rownames(res) = NULL
  res
}





