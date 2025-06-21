# Compute total log-likelihood with the forward algorithm
.loglik = function(.data, k1, a, Xchrom = FALSE, sex = NULL) {
  n = nrow(.data)
  logc = numeric(n)
  d = c(0, diff(.data$cm)) / 100

  # Initial state probabilities
  inits = c(nonIBD = 1 - k1, IBD = k1)

  emiss = as.matrix(.data[c("emission0","emission1")])

  # Forward - with logliks
  fwd = matrix(0, 2, n)
  alpha      = inits * emiss[1,]
  fwd[, 1]   = alpha / sum(alpha)
  logc[1]    = log(sum(alpha))

  for(i in 2:n) {
    tm = trans(d[i], k1, a)
    alpha = (tm %*% fwd[,i-1]) * emiss[i,]
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
#' @param prepped A logical indicating if the input data has undergone internal
#'   prepping. Can be ignored by most users.
#'
#' @returns A number: The total log-likelihood of the data under the HMM model.
#'
#' @examples
#' totalLoglik(cousinsDemo, k1 = 0.2, a = 5)
#'
#' @export
totalLoglik = function(data, ids = NULL, k1, a, prepped = FALSE) {

  .data = if(prepped) data else prepForHMM(data, ids = ids)

  Xchrom = attr(.data, "Xchrom")
  sex = attr(.data, "sex")

  chroms = unique.default(.data$chrom)
  logliks = numeric(length(chroms))

  # Loop through chroms
  for(i in seq_along(chroms)) {
    subdat = .data[.data$chrom == chroms[i], , drop  = FALSE]
    logliks[i] = .loglik(subdat, k1, a, Xchrom = Xchrom, sex = sex)
  }

  sum(logliks)
}
