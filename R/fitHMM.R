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
#' @param data Data frame with required columns `chrom`, `cm`, `a1` and `freq1`.
#' @param ids Genotype columns (default: last 2 columns).
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
#' @importFrom stats optim
#' @importFrom forrel ibdEstimate
#' @export
fitHMM = function(data, ids = NULL, thompson = TRUE, verbose = FALSE, ...) {

  .data = prepForHMM(data, ids = ids)

  # Initial parameter values
  k1_init = 0.5
  a_init = 50

  # Control parameters
  control = if(verbose) list(trace = 3, REPORT = 1) else control = list(...)

  # X chromosome
  Xchrom = all(.data$chrom == 23)
  if(Xchrom) {
    fnx = function(p) -totalLoglik(.data, k1 = p[1], a = p[2], Xchrom = TRUE, sex = getsex(.data))
    res = optim(c(k1 = k1_init, a = a_init), fnx, method = "L-BFGS-B",
                lower = c(0.001,0.001), upper = c(0.999, 100), control = control)
    return(as.list(res$par))
  }

  # Autosomal method 1: thompson optimisation of k1

  if(thompson) {
    s = asSingletons(.data, prep = FALSE)
    khat = ibdEstimate(s, verbose = FALSE)
    if(verbose)
      cat(sprintf("Initial kappa estimate: %.4g\n", khat))
    k1 = khat[[1, "k1"]]

    # Optimise `a` with k1 fixed
    fn1 = function(a) -totalLoglik(.data, k1 = k1, a = a)
    res = optim(a_init, fn1, method = "L-BFGS-B", lower = 0.001, upper = 100, control = control)
    return(list(k1 = k1, a = res$par))
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
