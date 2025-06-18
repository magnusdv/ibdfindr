library(pedsuite)
library(ibdsim2)
library(tibble)

x = cousinPed(1, removal = 1) |> setSNPs(snpData = FORCE)
ids = leaves(x)

sim = ibdsim(x, seed = 1)
findPattern(sim, list(carriers = ids))

y = profileSimIBD(x, sim, ids = ids)

g = getGenotypes(y, ids = ids)
stopifnot(identical(colnames(g), FORCE$MARKER))

cousinsDemo = tibble(FORCE, ID1 = g[1, ], ID2 = g[2, ]) |>
  add_column(CM = convertPos(chrom = FORCE$CHROM, Mb = FORCE$MB), .before = "A1")

usethis::use_data(cousinsDemo, overwrite = TRUE)


# Testing -----------------------------------------------------------------

pars = fitHMM(cousinsDemo)
post = ibdPosteriors(cousinsDemo, k1 = .2, a = pars$a)
segs = findSegments(cousinsDemo, k1 = .2, a = pars$a)
plotIBD(post, segs)
realisedKappa(sim, ids)
ibdEstimate(y, ids)
kappaIBD(x, ids)

