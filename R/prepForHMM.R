
prepForHMM = function(data, ids = NULL, err = 0, keepOld = FALSE) {

  # Check required columns
  lownames = tolower(names(data))
  required = c("chrom", "a1", "freq1")
  if(!all(required %in% lownames))
    stop2("Missing column: ", setdiff(required, lownames))

  if(!"cm" %in% lownames && !"mb" %in% lownames)
    stop2("Missing column: Either `cm` or `mb` must be present")

  # Make lower-case names for these
  names(data)[lownames == "chrom"] = "chrom"
  names(data)[lownames == "a1"] = "a1"
  names(data)[lownames == "freq1"] = "freq1"
  names(data)[lownames == "cm"] = "cm"
  names(data)[lownames == "mb"] = "mb"

  # Genotype columns (default: last 2 columns)
  if(is.null(ids))
    ids = utils::tail(names(data), 2)

  if(length(ids) != 2)
    stop2("Argument `ids` must have length exactly 2: ", ids)

  if(!all(ids %in% names(data)))
    stop2("Genotype column not found: ", setdiff(ids, names(data)))

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
  good = !is.na(data$chrom) & !is.na(data$cm %||% data$mb) & !is.na(data$a1) & !is.na(data$freq1)
  data = data[good, , drop = FALSE]

  # Convert to numeric chrom labels
  chromNum = chromNumber(data$chrom)
  badchroms = chromNum > 23 | chromNum < 1 | chromNum != as.integer(chromNum)
  if(any(badchroms))
    stop2("Illegal chromosomes: ", unique(data$chrom[badchroms]), "\nPlease remove these and try again")

  data$chrom = chromNum

  # Check for X chromosome
  Xchrom = all(chromNum == 23)
  if(!Xchrom && any(chromNum == 23))
    stop2("Cannot mix autosomal chromosomes and X chromosome. These must be processed separately.")

  # Convert mb to cm if needed
  if(!"cm" %in% names(data)) {
    data = mb2cm(data, chromCol = "chrom", mbCol = "mb", Xchrom = Xchrom)
  }

  # Deduce sex if X
  sex = if(Xchrom) getsex(data, ids = ids) else NULL

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
                              Xchrom = Xchrom, sex = sex, err = err)
  data$emission0 = emissionProbs[1,]
  data$emission1 = emissionProbs[2,]

  # Split in separate chroms
  data = split(data, data$chrom)

  # Return with attributes
  structure(data, ids = ids, Xchrom = Xchrom, sex = sex)
}

#' @importFrom ibdsim2 convertPos loadMap
mb2cm = function(data, chromCol, mbCol, Xchrom = FALSE) {

  if(!Xchrom)
    cm = ibdsim2::convertPos(chrom = data[[chromCol]], Mb = data[[mbCol]], sex = "average")
  else {
    mapX = ibdsim2::loadMap(chrom = "X")[[1]]
    cm = ibdsim2::convertPos(Mb = data[[mbCol]], map = mapX, sex = "female")
  }

  data$cm = cm
  data
}
