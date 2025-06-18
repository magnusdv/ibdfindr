stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Normalise vector to sum 1
nrm = function(x) x/sum(x)

# Return IBS states of two genotype vectors
ibsState = function(g1, g2) {
  ibs = rep(1L, length(g1))
  ibs[g1 == g2] = 2L
  ibs[g1 == "1/1" & g2 == "2/2"] = 0L
  ibs[g1 == "2/2" & g2 == "1/1"] = 0L
  ibs
}


getsex = function(data, ids = c("g1", "g2")) {
  sapply(ids, function(id) if(all(nchar(data[[id]]) == 1)) 1 else 2)
}

# Convert genotypes to numeric, e.g. "AC" -> "1/2"
g2num = function(g, allele1, Xchrom = FALSE) {
  # Simple return if X and male
  if(Xchrom && all(nchar(g) == 1))
    return(ifelse(g == allele1, "1", "2"))

  # Otherwise
  hasSep = nchar(g[1]) == 3
  a = ifelse(substr(g, 1, 1) == allele1, 1, 2)
  if(hasSep)
    b = ifelse(substr(g, 3, 3) == allele1, 1, 2)
  else
    b = ifelse(substr(g, 2, 2) == allele1, 1, 2)

  gnum = paste(a, b, sep = "/")
  gnum[gnum== "2/1"] = "1/2"
  gnum
}

# Matrix with probabilities of observed emissions in the given dataset
emissionMat = function(freq1, g1, g2, Xchrom = FALSE, sex = NULL) {
  em = vapply(seq_along(freq1), function(i) {
    emat = emission(freq1[i], Xchrom = Xchrom, sex = sex)
    emat[g1[i], g2[i], ]
  }, FUN.VALUE = numeric(2))

  colnames(em) = names(freq1)
  em
}

# Numeric chromosome labels (for sorting)
chromNumber = function(x) {
  # Strip prefixes
  labs = sub("^(chrom|chr)", "", x, ignore.case = TRUE)

  num  = suppressWarnings(as.numeric(labs))
  num[labs %in% c("X","x")] = 23
  num[labs %in% c("Y","y")] = 24

  # Remaining NAs go after, in original relative order
  na_idx = is.na(num)
  num[na_idx] = max(num, na.rm = TRUE) + seq_len(sum(na_idx))

  num
}

#' @importFrom pedtools singletons setSNPs setGenotype
asSingletons = function(data, ids = NULL, prep = TRUE) {
  ids = ids %||% utils::tail(names(data), 2)

  if(prep)
    data = prepForHMM(data, ids = ids)

  snpData = data.frame(CHROM = data$chrom, MARKER = data$marker, MB = data$cm,
                       A1 = "1", A2 = "2", FREQ1 = data$freq1)
  singletons(ids) |>
    setSNPs(snpData = snpData) |>
    setGenotype(ids = ids[1], geno = data$g1) |>
    setGenotype(ids = ids[2], geno = data$g2)
}
