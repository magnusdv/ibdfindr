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

  sex = if(Xchrom) getsex(data) else NULL

  # Sort
  data = data[order(data$chrom, data$cm), , drop = FALSE]

  # Numeric genotype columns
  data$g1 = g2num(data[[ids[1]]], data$a1, Xchrom = Xchrom)
  data$g2 = g2num(data[[ids[2]]], data$a1, Xchrom = Xchrom)

  if(!keepOld)
    data = data[, c("chrom", "marker", "cm", "freq1", "g1", "g2"), drop = FALSE]

  # IBS status
  data$ibs = ibsState(data$g1, data$g2, Xchrom = Xchrom)

  # Emission probs conditioned on IBD 0 and 1
  emissionProbs = emissionMat(data$freq1, g1 = data$g1, g2 = data$g2,
                              Xchrom = Xchrom, sex = sex)
  data$emission0 = emissionProbs[1,]
  data$emission1 = emissionProbs[2,]

  # Return with attributes
  structure(data, ids = ids, Xchrom = Xchrom, sex = sex)
}
