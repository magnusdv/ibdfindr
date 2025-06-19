library(pedsuite)
library(ibdsim2)
library(tibble)

# Pedigree
x = nuclearPed(2) |> setSNPs(snpData = XFORCE)
ids = 3:4

# Simulate IBD pattern
mapX = loadMap(chrom = 23)[[1]]
sim = ibdsim(x, seed = 1, map = mapX, verbose = F)

# True IBD sharing
trueSegs = findPattern(sim, list(carriers = ids))
trueSegs
haploDraw(x, sim)

# Simulate genotypes
y = profileSimIBD(x, sim, ids = ids, seed = 1729)
g = getGenotypes(y, ids = ids)

# Dataset with annotation and genotypes
brothersX = tibble(XFORCE, ID1 = g[1, ], ID2 = g[2, ])

# Add centiMorgan positions
cm = convertPos(Mb = XFORCE$MB, map = mapX, sex = "female")
brothersX = brothersX |> add_column(CM = cm, .before = "A1")

# Save dataset
usethis::use_data(brothersX, overwrite = TRUE)


# Testing -----------------------------------------------------------------

library(ibdfindr)
pars = fitHMM(brothersX)
pars
post = ibdPosteriors(brothersX, k1 = pars$k1, a = pars$a)
segs = findSegments(brothersX, k1 = pars$k1, a = pars$a)
plotIBD(post, segs)

# With true segments
plotIBD(post, segs) +
  geom_segment(data = trueSegs, col = "blue", linewidth = 1.5,
               aes(x = startCM, xend = endCM, y = .2, yend = .2))
