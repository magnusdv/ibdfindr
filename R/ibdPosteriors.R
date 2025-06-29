# Compute posterior IBD probabilities for a single chromosome
.ibdPosteriors = function(.data, k1, a, Xchrom = FALSE, sex = NULL) {
  if(length(unique.default(.data$chrom)) > 1)
    stop("Function `.ibdPosterior()` expects a single chromosome.")

  n = nrow(.data)
  d = c(0, diff(.data$cm)) / 100

  trans = transitionProbs(d, k1, a)
  emiss = as.matrix(.data[c("emission0","emission1")])

  # Forward-backward algorithm ----------------------------------------------

  # Initial state probabilities
  inits = c(nonIBD = 1 - k1, IBD = k1)

  # Forward (with logliks)
  fwd = matrix(0, 2, n)
  fwd[, 1] = (inits * emiss[1,]) |> nrm()

  idxF = if(n>1) 2:n else integer(0)
  for(i in idxF) {
    tt = trans[i, ] * fwd[, i-1]
    fwd[,i] = (c(tt[1] + tt[2], tt[3] + tt[4]) * emiss[i,]) |> nrm()
  }

  # Backward
  bwd = matrix(1, 2, n)
  idxB = if(n>1) (n-1):1 else integer(0)
  for(i in idxB) {
    tt = trans[i+1, ] * emiss[i+1,] * bwd[,i+1]
    bwd[,i] = c(tt[1] + tt[2], tt[3] + tt[4]) |> nrm()
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
#' @param prepped A logical indicating if the input data has undergone internal
#'   prepping. Can be ignored by most users.
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
ibdPosteriors = function(data, ids = NULL, k1, a, prepped = FALSE, verbose = FALSE) {

  cat("Calculating IBD posteriors...\n")

  .data = if(prepped) data else prepForHMM(data, ids = ids)

  Xchrom = attr(.data, "Xchrom")
  sex = attr(.data, "sex")

  # Loop through chroms
  dflist = lapply(.data, function(chrdat) {
    post = .ibdPosteriors(chrdat, k1, a, Xchrom = Xchrom, sex = sex)
    chrdat$post = post
    chrdat
  })

  res = do.call(rbind, dflist)
  rownames(res) = NULL
  res
}





