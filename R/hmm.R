
# -------------------------------------------------------------------------
# Elements of HMM for IBD status ------------------------------------------
# -------------------------------------------------------------------------


# Transition probabilities for two-state CT-HMM
# (similar to Leutenegger et al 2003)
# Output: Matrix with one row per step
# Columns: 0->0, 0->1, 1->0, 1->1, where 0 = nonIBD, 1 = IBD
transitionProbs = function(dvec, k1, a) {
  A = exp(-a * dvec)
  cbind((1-A) * (1 - k1) + A,
        (1-A) * k1,
        (1-A) * (1 - k1),
        (1-A) * k1 + A)
}


# Emission probs P(G1, G2 | IBD) for given allele freqs (p, q).
# Returns array with IBD = 0/1 in third dimension.
emission = function(p, Xchrom = FALSE, sex = NULL, err = 0) {
  fr = c(p, q <- 1-p)
  hw = c(p^2, 2*p*q, q^2)
  a = c("1", "2")
  g = c("1/1", "1/2", "2/2")
  dim3 = c("nonIBD", "IBD")

  # Autosomal or X female/female
  if(!Xchrom || all(sex == 2)) {
    mat0 = hw %o% hw
    mat1 = matrix(c(p^3, p^2*q, 0,
                    p^2*q, p*q, p*q^2,
                    0, p*q^2, q^3), 3, 3)
    E = array(c(mat0, mat1), dim = c(3, 3, 2), dimnames = list(g, g, dim3))
  }
  else if(all(sex == 1)) { # X: male/male
    mat0 = fr %o% fr
    mat1 = diag(fr)
    E = array(c(mat0, mat1), dim = c(2, 2, 2), dimnames = list(a, a, dim3))
  }
  else if(sum(sex) == 3) { # X: female/male or vice versa)
    mat0 = hw %o% fr
    mat1 = matrix(c(p^2, p*q, 0, 0, p*q, q^2), 3, 2)
    if(sex[1] == 2)
      E = array(c(mat0, mat1), dim = c(3, 2, 2), dimnames = list(g, a, dim3))
    else
      E = array(c(t(mat0), t(mat1)), dim = c(2, 3, 2), dimnames = list(a, g, dim3))
  }
  else
    stop2("Illegal `sex`:", sex)

  if(err > 0)
    E[,,"IBD"] = (1-err) * E[,,"IBD"] + err * E[,,"nonIBD"]

  E
}


# Matrix with probabilities of observed emissions in the given dataset
emissionMat = function(freq1, g1, g2, Xchrom = FALSE, sex = NULL, err = 0) {
  if(length(err) != 1 || err < 0 || err > 1)
    stop2("`err` must be a single number between 0 and 1")

  em = vapply(seq_along(freq1), function(i) {
    emat = emission(freq1[i], Xchrom = Xchrom, sex = sex, err = err)
    emat[g1[i], g2[i], ]
  }, FUN.VALUE = numeric(2))

  colnames(em) = names(freq1)
  em
}
