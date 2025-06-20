#' Fit a Hidden Markov Model to genotype data
#'
#' This function fits a continuous-time HMM to the given genotype data. It does
#' so by optimising the parameters `k1` (the probability of being in an IBD
#' state) and the transition rate `a` to maximise the total log-likelihood.
#'
#' If `thompson = TRUE` (default), the parameter `k1` is estimated first, using
#' the standard maximum-likelihood approach for estimating pairwise relatedness
#' coefficients (first described by Thompson 1975). The actual work is done by a
#' call to [forrel::ibdEstimate()]. The parameter `a` is subsequently estimated
#' with `optim`.
#'
#' If `thompson = FALSE`, both parameters `k1` and `a` are optimised together.
#' (This is the only method for X chromosomal data.)
#'
#' @param data Data frame with required columns `chrom`, `cm`, `a1` and `freq1`
#'   (case insensitive).
#' @param ids Genotype columns (default: last 2 columns).
#' @param k1,a Numeric HMM parameters. Supplying a value fixes the parameter; if
#'   NULL (default), the parameter is estimated.
#' @param thompson A logical indicating the optimisation method. (See Details.)
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
#' @importFrom stats optim optimise
#' @importFrom forrel ibdEstimate
#' @export
fitHMM = function(data, ids = NULL, k1 = NULL, a = NULL, thompson = TRUE, verbose = FALSE, ...) {

  .data = prepForHMM(data, ids = ids)

  # Initial parameter values
  k1_init = k1 %||% 0.5
  a_init = a %||% 5

  # X and sex
  Xchrom = all(.data$chrom == 23)
  sex = if(Xchrom) getsex(.data) else NULL

  if(verbose) {
    if(Xchrom)
      cat(sprintf("Chromosome type: X (%s)\n", paste(c("male", "female")[sex], collapse="/")))
    else
      cat("Chromosome type: autosomal\n")
  }
  # If `k1` or `a` are provided: estimate the other -------------------------

  if(!is.null(k1) && !is.null(a)) {
    if(verbose) cat("Nothing to do; parameters provided: k1 =", k1, ", a =", a, "\n")
    return(list(k1 = k1, a = a))
  }
  if(!is.null(k1) && is.null(a)) {
    if(verbose) cat("Optimising `a` conditional on k1 =", k1, "\n")
    fn01 = function(a) -totalLoglik(.data, k1 = k1, a = a, Xchrom = Xchrom, sex = sex)
    res = optimise(fn01, interval = c(0.001, 100), tol = 0.001)
    return(list(k1 = k1, a = res$minimum))
  }
  if(is.null(k1) && !is.null(a)) {
    if(verbose) cat("Optimising `k1` conditional on a =", a, "\n")
    fn02 = function(k1) -totalLoglik(.data, k1 = k1, a = a, Xchrom = Xchrom, sex = sex)
    res = optimise(fn02, interval = c(0.001, 0.999), tol = 0.001)
    return(list(k1 = res$minimum, a = a))
  }


  # Otherwise: Estimate both ------------------------------------------------

  # Control parameters for `optim`
  control = list(...)

  # X chromosome
  if(Xchrom) {
    fnx = function(p) -totalLoglik(.data, k1 = p[1], a = p[2], Xchrom = TRUE, sex = sex)
    res = optim(c(k1 = k1_init, a = a_init), fnx, method = "L-BFGS-B",
                lower = c(0.001,0.001), upper = c(0.999, 100), control = control)
    return(as.list(res$par))
  }

  if(verbose) {
    cat("Thompson estimation:", thompson, "\n")
  }
  # Autosomal method 1: thompson optimisation of k1

  if(thompson) {
    s = asSingletons(.data, prep = FALSE)
    khat = ibdEstimate(s, verbose = FALSE)[1, 4:6] |> as.numeric()
    k1 = khat[2]
    if(verbose) {
      cat(sprintf("Estimated kappa: (k0, k1, k2) = (%s)\n", toString(round(khat,3))))
      cat("Optimising `a` conditional on k1 =", round(khat[2],3), "\n")
    }

    # Optimise `a` with k1 fixed
    fn1 = function(a) -totalLoglik(.data, k1 = k1, a = a)
    res = optimise(fn1, interval = c(.01, 100), tol = 0.001)  # a bit faster than `optim`
    return(list(k1 = k1, a = res$minimum))
  }

  # Autosomal method 2: optimise k1 and a together

  fn2 = function(p) -totalLoglik(.data, k1 = p[1], a = p[2])
  res = optim(c(k1 = k1_init, a = a_init), fn2, method =  "L-BFGS-B",
              lower = c(0.001,0.001), upper = c(0.999, 100), control = control)
  as.list(res$par)
}



### Attempt at Baum-Welch. Not working properly.
### Not vital to get this working; thompson-optimisation as above probably better (!?)
fitHMM_BW = function(data, ids = NULL, maxIter = 20, tol = 1e-6, verbose = T) {

  data = prepForHMM(data, ids = ids)

  # Init parameters
  k1 = 0.2
  a = 10

  if(verbose) cat(sprintf("Start: k1=%.4g  a=%.4g  logLik=%.3f\n",
                          k1, a, totalLoglik(data, k1 = k1, a = a)))

  ntot = nrow(data)
  chromIdx = split(seq_len(ntot), data$chrom)
  p = data$freq1
  finLog = \(x) log(pmax(x, .Machine$double.eps))
  aRange = c(0.01, 100)

  for(it in 1:maxIter) {
    gSum = 0; nTot = 0; logL.old = 0
    xiAll = list(); dAll = numeric()

    ## ----- E-step across chromosomes -----
    for(idx in chromIdx) {
      n  = length(idx)
      d  = c(0, diff(data$cm[idx]))/100
      em = sapply(idx, \(i) emission(p[i])[data$g1[i], data$g2[i], ])

      # Forward
      fwd = bwd = matrix(0, 2, n)
      logc = numeric(n)
      alpha = c(1-k1, k1) * em[, 1]
      logc[1] = log(sum(alpha))
      fwd[, 1] = alpha/sum(alpha)
      for(i in 2:n) {
        alpha = (trans(d[i], k1, a) %*% fwd[, i-1]) * em[,i]
        logc[i] = log(sum(alpha))
        fwd[, i] = alpha/sum(alpha)
      }

      bwd[, n] = 1
      for(i in (n-1):1)
        bwd[, i] = (trans(d[i+1], k1, a) %*% (em[,i+1] * bwd[,i+1])) |> nrm()

      # Posteriors
      gamma = fwd * bwd
      gamma = gamma / rep(colSums(gamma), each = 2)


      xi = array(0, c(2, 2, n-1))
      for(i in 2:n) {
        tm = trans(d[i], k1, a)
        xi[,,i-1] = ((fwd[, i-1] %o% (em[, i] * bwd[, i])) * tm) |> nrm()
      }
      gSum = gSum + sum(gamma[2, ])
      nTot = nTot + n
      xiAll = c(xiAll, asplit(xi, 3))
      dAll = c(dAll, d[-1])
      logL.old = logL.old + sum(logc)
    }

    ## ----- M-step -----
    k1.new = gSum/nTot
    ll.a = function(aa)
      -sum(vapply(seq_along(dAll), \(j) sum(xiAll[[j]] * finLog(trans(dAll[j], k1.new, aa))), 0.0))
    a.new = stats::optimize(ll.a, aRange)$minimum

    ## likelihood under updated params
    logL.new = totalLoglik(data, k1 = k1.new, a = a.new)
    if(verbose)
      cat(sprintf("it=%d  k1=%.4g  a=%.4g  logLik(old)=%.3f  logLik(new)=%.3f\n",
                  it, k1.new, a.new, logL.old, logL.new))

    if(abs(logL.new - logL.old) < tol) {
      k1 = k1.new
      a = a.new
      break
    }
    k1 = k1.new
    a = a.new
  }

  list(k1 = k1, a = a, logLik = totalLoglik(data, k1 = k1, a = a), iter = it)
}
