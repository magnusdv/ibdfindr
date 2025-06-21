
# -------------------------------------------------------------------------
# Elements of HMM for IBD status ------------------------------------------
# -------------------------------------------------------------------------


# Transition matrix for two-state CT-HMM
trans = function(d, k1, a) {
  e = exp(-a * d)
  matrix(c((1 - e)*(1 - k1) + e, (1 - e)*k1,
           (1 - e)*(1 - k1),     (1 - e)*k1 + e),
         2, 2, byrow = TRUE, dimnames = rep(list(c("nonIBD", "IBD")), 2))
}


# Emission probs P(G1, G2 | IBD) for given allele freqs (p, q).
# Returns array with IBD = 0/1 in third dimension.
# TODO: Possible to simplify using IBS status?
emission = function(p, Xchrom = FALSE, sex = NULL) {
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
    array(c(mat0, mat1), dim = c(3, 3, 2), dimnames = list(g, g, dim3))
  }
  else if(all(sex == 1)) { # X: male/male
    mat0 = fr %o% fr
    mat1 = diag(fr)
    array(c(mat0, mat1), dim = c(2, 2, 2), dimnames = list(a, a, dim3))
  }
  else if(sum(sex) == 3) { # X: female/male or vice versa)
    mat0 = hw %o% fr
    mat1 = matrix(c(p^2, p*q, 0, 0, p*q, q^2), 3, 2)
    if(sex[1] == 2)
      array(c(mat0, mat1), dim = c(3, 2, 2), dimnames = list(g, a, dim3))
    else
      array(c(t(mat0), t(mat1)), dim = c(2, 3, 2), dimnames = list(a, g, dim3))
  }
  else
    stop2("Illegal `sex`:", sex)
}


