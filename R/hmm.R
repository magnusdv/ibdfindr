
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


prepForHMM = function(data, ids = NULL, keepOld = FALSE) {

  # Check required columns
  required = c("chrom", "cm", "a1", "freq1")
  lownames = tolower(names(data))
  if(!all(required %in% lownames))
    stop2("Missing column: ", setdiff(required, lownames))

  # Use lower-case names for these
  names(data)[match(required, lownames)] = required

  # Genotype columns (default: last 2 columns)
  if(is.null(ids))
    ids = utils::tail(names(data), 2)

  if(length(ids) != 2)
    stop2("Expected 2 genotype columns, got ", length(ids))

  if(!all(ids %in% names(data)))
    stop2("Genotype columns not found: ", setdiff(ids, names(data)))

  if(any(tolower(ids) %in% required))
    stop2("Illegal genotype column: ", intersect(tolower(ids), required))

  # Marker column
  midx = match(c("marker","snp","name","rs"), lownames, nomatch = 0L)
  midx = midx[midx > 0]
  if(length(midx))
    names(data)[midx[1]] = "marker"
  else
    data$marker = paste0("m", seq_len(nrow(data)))

  names(data$freq1) = data$marker

  # Remove missing data
  data = data[!is.na(data$chrom) & !is.na(data$cm) & !is.na(data$a1) & !is.na(data$freq1), , drop = FALSE]

  # Convert to numeric chrom labels
  chromNum = chromNumber(data$chrom)
  badchroms = chromNum > 23 | chromNum < 1 | chromNum != as.integer(chromNum)
  if(any(badchroms))
    stop2("Illegal chromosomes: ", unique(data$chrom[badchroms]), "\nPlease remove these and try again")

  data$chrom = chromNum

  # Check for X chromosome
  Xchrom = all(chromNum == 23)
  if(!Xchrom && any(chromNum == 23))
    stop2("Cannot mix chromosome 23 (X) with autosomes")

  # Sort
  data = data[order(data$chrom, data$cm), , drop = FALSE]

  data$g1 = g2num(data[[ids[1]]], data$a1, Xchrom = Xchrom)
  data$g2 = g2num(data[[ids[2]]], data$a1, Xchrom = Xchrom)
  data$ibs = ibsState(data$g1, data$g2)

  if(!keepOld)
    data = data[, c("chrom", "marker", "cm", "freq1", "g1", "g2", "ibs"), drop = FALSE]

  data
}
