#' Fit a Hidden Markov Model to genotype data
#'
#' This function fits a continuous-time HMM to the provided genotype data, by
#' optimising the parameters `k1` (the probability of being in an IBD state) and
#' `a` (the transition rate) to maximise the total log-likelihood.
#'
#' By default (`thompson = FALSE`) both parameters `k1` and `a` are optimised
#' together, using [stats::optimise()].
#'
#' If `thompson = TRUE`, then `k1` is estimated first, using the
#' maximum-likelihood approach for pairwise relatedness coefficients described
#' by Thompson (1975). (Note that, although this method was originally developed
#' for unlinked markers, it yields unbiased estimates also with linked markers.)
#' The estimation of `k1` is performed internally by calling
#' [forrel::ibdEstimate()]. Subsequently, the parameter `a` is estimated
#' conditional on the `k1` value.
#'
#' @param data Data frame with required columns `chrom`, `cm`, `a1` and `freq1`
#'   (case insensitive).
#' @param ids Genotype columns (default: last 2 columns).
#' @param k1,a Numeric HMM parameters. Supplying a value fixes the parameter; if
#'   NULL (default), the parameter is estimated.
#' @param err Error rate; a single number in `[0,1]` (default: 0).
#' @param method A character string indicating the optimisation method to use.
#' @param thompson A logical indicating the optimisation method. (See Details.)
#' @param prepped A logical indicating if the input data has undergone internal
#'   prepping. Can be ignored by most users.
#' @param verbose A logical indicating whether to print information during the
#'   optimisation.
#' @param ... Additional arguments passed to the `control` parameter of
#'   [stats::optim()].
#'
#' @returns A list containing the fitted parameters `k1` and `a`, and some
#'   additional information about the optimisation.
#'
#' @seealso [totalLoglik()], [forrel::ibdEstimate()]
#'
#' @examples
#' \dontrun{
#' fitHMM(cousinsDemo)
#' }
#'
#' @importFrom stats optim optimise plogis qlogis
#' @importFrom forrel ibdEstimate
#' @export
fitHMM = function(data, ids = NULL, k1 = NULL, a = NULL, err = 0,
                  method = "L-BFGS-B", thompson = FALSE, prepped = FALSE,
                  verbose = FALSE, ...) {

  .data = if(prepped) data else prepForHMM(data, ids = ids, err = err)

  # Parameter bounds
  k1_min = 0.01
  k1_max = 0.99
  a_min = 0.01
  a_max = 30

  # Initial parameter values
  k1_init = k1 %||% 0.2
  a_init = a %||% 5

  # X and sex
  Xchrom = attr(.data, "Xchrom")
  sex = attr(.data, "sex")

  if(verbose) {
    if(Xchrom)
      chrtype = sprintf("X (%s)", paste(c("male", "female")[sex], collapse="/"))
     else
       chrtype = "autosomal"

    cat("Chromosome type:", chrtype, "\n")
    cat("Fitting HMM parameters...\n")
  }

  # Optional: Thompson estimation of k1
  if(thompson && is.null(k1)) {

    if(Xchrom)
      stop2("Thompson estimation of `k1` is not implemented for X chromosome data.")

    s = asSingletons(.data, prepped = TRUE)
    khat = ibdEstimate(s, verbose = FALSE)[1, 4:6] |> as.numeric()
    k1 = khat[2]

    if(verbose)
      cat(sprintf("  Thompson estimate: (k0, k1, k2) = (%s)\n",
                  toString(round(khat,3))))
  }

  if(!is.null(k1) && !pedtools:::isNumber(k1, 0, 1))
    stop2("`k1` must be a number in the interval `[0, 1)`: ", k1)

  if(!is.null(a) && !pedtools:::isNumber(a, 0))
      stop2("`a` must be a positive number: ", a)

  # Estimate jointly if neither are provided

  if(is.null(k1) && is.null(a)) {
    meth = method %||% "L-BFGS-B"
    if(verbose) cat("  Optimising `k1` and `a` jointly; method:", meth, "\n")

    if(meth == "L-BFGS-B") {
      fn1 = function(p) -totalLoglik(.data, k1 = p[1], a = p[2], prepped = TRUE)
      res = optim(c(k1 = k1_init, a = a_init), fn1, method =  "L-BFGS-B",
                lower = c(k1_min, a_min), upper = c(k1_max, a_max), control = list(...))
      respar = as.list(res$par); #print(res$count)
      k1 = respar$k1
      a = respar$a
    }
    else {
      k1_init = plogis(k1_init)
      a_init = log(a_init)
      fn1b = function(q)
        -totalLoglik(.data, k1 = plogis(q[1]), a = exp(q[2]), prepped = TRUE)
      res = optim(c(k1 = k1_init, a = a_init), fn1b, method = method, control = list(...))
      respar = as.list(res$par); #print(res$count)
      k1 = plogis(respar$k1)
      a = exp(respar$a)
    }
  }
  else if(!is.null(k1) && is.null(a)) {
    if(verbose) cat("  Optimising `a` conditional on k1 =", k1, "\n")

    fn2 = function(a) -totalLoglik(.data, k1 = k1, a = a, prepped = TRUE)
    res = optimise(fn2, interval = c(a_min, a_max), tol = 0.001)
    a = res$minimum
  }
  else if(is.null(k1) && !is.null(a)) {
    if(verbose) cat("  Optimising `k1` conditional on a =", a, "\n")

    fn3 = function(k1) -totalLoglik(.data, k1 = k1, a = a, prepped = TRUE)
    res = optimise(fn3, interval = c(k1_max, k1_min), tol = 0.001)
    k1 = res$minimum
  }

  # Report and return
  if(verbose) cat(sprintf("  k1 = %.3f, a = %.3f\n", k1, a))
  list(k1 = k1, a = a)
}


# Estimate HMM parameters (k1, a) from pedigree
paramsFromPed = function(ped, ids = leaves(ped)) {

  # k1: The pedigree coefficient kappa_1
  inb = ribd::inbreeding(ped, ids = ids)
  if(any(inb > 0))
    stop2("Some of the indicated individuals are inbred; cannot proceed")
  kap = ribd::kappaIBD(ped, ids = ids, inbredAction = 0)
  if(kap[3] > 0)
    stop2("Pedigree relationship is not unilineal; cannot proceed")

  k1 = kap[2]

  # a: Estimate via two-locus k_11 over a range of recombination ratios
  cm = 0:20
  rho = 0.5 * (1 - exp(-cm/50))

  k11 = unname(sapply(rho, function(r)
    ribd::twoLocusIBD(ped, ids, rho = r, coefs = "k11")))

  # Non-linear fit for a:
  nlf = function(a) {
    pred = k1^2 + (1 - k1)*k1*exp(-a * cm/100)
    sum((k11 - pred)^2)
  }

  fit = stats::optimize(nlf, c(0, 20), tol = 1e-3, maximum = FALSE)
  a = round(fit$minimum, 3)

  # Return parameters
  list(k1 = k1, a = a)
}
