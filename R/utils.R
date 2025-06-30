stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Normalise vector to sum 1
nrm = function(x) x/sum(x)

# Return IBS states of two genotype vectors
ibsState = function(g1, g2, Xchrom = FALSE) {
  ibs = rep(1L, length(g1))

  if(!Xchrom) {
    ibs[g1 == g2] = 2L
    ibs[g1 == "1/1" & g2 == "2/2"] = 0L
    ibs[g1 == "2/2" & g2 == "1/1"] = 0L
  }
  else { # TODO: not very clever
    male1 = all(nchar(g1) == 1)
    male2 = all(nchar(g2) == 1)
    if(male1 && male2) {
      ibs[g1 != g2] = 0L
    }
    else if(male1 && !male2) {
      ibs[g1 == "1" & g2 == "2/2"] = 0L
      ibs[g1 == "2" & g2 == "1/1"] = 0L
    }
    else if(!male1 && male2) {
      ibs[g1 == "1/1" & g2 == "2"] = 0L
      ibs[g1 == "2/2" & g2 == "1"] = 0L
    }
    else {
      ibs[g1 == g2] = 2L
      ibs[g1 == "1/1" & g2 == "2/2"] = 0L
      ibs[g1 == "2/2" & g2 == "1/1"] = 0L
    }
  }
  ibs
}

# Deduce sex from X genotypes: Hemizygous -> male
getsex = function(data, ids = c("g1", "g2")) {
  if(!all(ids %in% names(data)))
    stop2("Genotype column not found: ", setdiff(ids, names(data)))

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


asSingletons = function(data, ids = NULL, prepped = FALSE) {

  .data = if(prepped) data else prepForHMM(data, ids = ids, keepOld = TRUE)

  ids = attr(.data, "ids")
  sex = attr(.data, "sex")

  .data = do.call(rbind, .data)

  mb = if("MB" %in% names(.data)) .data$MB else .data$cm
  gcols = .data[paste0("g", seq_along(ids))]
  names(gcols) = ids

  snpData = data.frame(CHROM = .data$chrom,
                       MARKER = .data$marker,
                       MB = mb,
                       A1 = "1", A2 = "2",
                       FREQ1 = .data$freq1,
                       gcols,
                       check.names = FALSE)

  singletons(ids, sex = sex %||% 1) |>
    setSNPs(snpData = snpData)
}


#' @importFrom ibdsim2 loadMap convertPos
getSNPdata = function(x, ids = NULL) {
  snpmap = getMap(x)

  chrs = unique.default(snpmap$CHROM)
  Xchrom = length(chrs) == 1 && chrs %in% c("X", 23)
  gmap = loadMap("decode19", chrom = chrs)

  snpmap$CM = convertPos(chrom = snpmap$CHROM, Mb = snpmap$MB, map = gmap,
                         sex = if(Xchrom) "female" else "average")

  als = getLocusAttributes(x, attribs = "alleles", simplify = TRUE)
  if(!length(als))
    stop2("No attached markers")
  if(!all(lengths(als) == 2))
    stop2("Some markers are not SNPs")

  alsmat = do.call(rbind, als)
  colnames(alsmat) = c("A1", "A2")

  freq1 = sapply(x$MARKERS, function(m) attr(m, "afreq")[1])

  ids = ids %||% typedMembers(x)
  g = getGenotypes(x, ids = ids) |> t.default()

  cbind(snpmap, alsmat, FREQ1 = freq1, g)
}
