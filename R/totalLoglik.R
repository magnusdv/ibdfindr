# Compute total log-likelihood with the forward algorithm
.loglik = function(.data, k1, a, Xchrom = FALSE, sex = NULL) {
  n = nrow(.data)
  logc = numeric(n)
  d = c(0, diff(.data$cm)) / 100

  # Initial state probabilities
  inits = c(nonIBD = 1 - k1, IBD = k1)

  trans = transitionProbs(d, k1, a)
  emiss = cbind(.data$emission0, .data$emission1)

  # Forward - with logliks
  fwd = matrix(0, 2, n)
  alpha = inits * emiss[1,]
  fwd[, 1] = alpha / sum(alpha)
  logc[1] = log(sum(alpha))

  if(n == 1)
    return(logc[1])

  for(i in 2:n) {
    tt = trans[i, ] * fwd[, i-1]
    alpha = c(tt[1] + tt[2], tt[3] + tt[4]) * emiss[i,]
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
#' @param err Error rate; a single number in `[0,1]` (default: 0).
#' @param prepped A logical indicating if the input data has been internally
#'   processed. Can be ignored by most users.
#'
#' @returns A number: The total log-likelihood of the data under the HMM model.
#'
#' @examples
#' totalLoglik(cousinsDemo, k1 = 0.2, a = 5)
#'
#' @export
totalLoglik = function(data, ids = NULL, k1, a, err = 0, prepped = FALSE) {

  .data = if(prepped) data else prepForHMM(data, ids = ids, err = err)

  Xchrom = attr(.data, "Xchrom")
  sex = attr(.data, "sex")

  logliks = lapply(.data, function(chrdat)
    .loglik(chrdat, k1, a, Xchrom = Xchrom, sex = sex)
  )

  # Attempt with purrr::in_parallel()
  # logliks = .data |> map(in_parallel(
  #   \(chrdat) .loglik(chrdat, k1, a, Xchrom = Xchrom, sex = sex),
  #   k1 = k1, a = a, Xchrom = Xchrom, sex = sex, .loglik = ibdfindr:::.loglik))

  sum(unlist(logliks))
}
