# Compute total log-likelihood with the forward algorithm
.loglik = function(.data, k1, a, Xchrom = FALSE, sex = NULL) {
  n = nrow(.data)
  logc = numeric(n)
  d = c(0, diff(.data$cm)) / 100

  # Matrix (2*n) of emission probabilities
  emProb = emissionMat(.data$freq1, g1 = .data$g1, g2 = .data$g2, Xchrom = Xchrom, sex = sex)

  # Initial state probabilities
  inits = c(nonIBD = 1 - k1, IBD = k1)

  # Forward - with logliks
  fwd = matrix(0, 2, n)
  alpha      = inits * emProb[, 1]
  fwd[, 1]   = alpha / sum(alpha)
  logc[1]    = log(sum(alpha))

  for(i in 2:n) {
    tm = trans(d[i], k1, a)
    alpha = (tm %*% fwd[,i-1]) * emProb[,i]
    logc[i] = log(sum(alpha))
    fwd[,i] = alpha / sum(alpha)
  }

  # Return total loglik
  sum(logc)
}


#' Total log-likelihood for observed data
#'
#' This function computes the total log-likelihood of the observed data, under
#' the hidden Markov model. It is mainly for internal use, especially
#' [fitHMM()].
#'
#' @param data Data frame with required columns `chrom`, `cm`, `a1` and `freq1`.
#' @param ids Genotype columns (ignored unless `prep = TRUE`).
#' @param k1,a HMM parameters.
#' @param Xchrom A logical, by default FALSE.
#' @param sex A numeric of length 2, with entries 1 (=male) or 2 (=female).
#'   Ignored unless `Xhcrom = TRUE`.
#'
#' @returns A number: The total log-likelihood of the data under the HMM model.
#'
#' @examples
#' totalLoglik(cousinsDemo, k1 = 0.2, a = 5)
#'
#' @export
totalLoglik = function(data, ids = NULL, k1, a, Xchrom = NULL, sex = NULL) {
  prep = !all(c("chrom", "cm", "freq1", "g1", "g2") %in% names(data))
  if(prep)
    data = prepForHMM(data, ids = ids)

  chroms = unique.default(data$chrom)
  Xchrom = Xchrom %||% length(chroms) == 1 && chroms == 23
  if(Xchrom && is.null(sex))
    stop2("Sex info missing")

  logliks = numeric(length(chroms))

  # Loop through chroms
  for(i in seq_along(chroms)) {
    subdat = data[data$chrom == chroms[i], , drop  = FALSE]
    logliks[i] = .loglik(subdat, k1, a, Xchrom = Xchrom, sex = sex)
  }

  sum(logliks)
}
