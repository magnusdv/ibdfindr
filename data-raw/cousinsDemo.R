library(pedsuite)
library(ibdsim2)
library(tibble)

# Pedigree
x = cousinPed(1, removal = 1) |> setSNPs(snpData = FORCE)
ids = leaves(x)

# Simulate IBD pattern
sim = ibdsim(x, seed = 1)

# True IBD sharing
trueSegs = findPattern(sim, list(carriers = ids))
karyoHaploid(trueSegs)

# Simulate genotypes
y = profileSimIBD(x, sim, ids = ids, seed = 1729)

# Better coding of dataset
#cousinsDemo = getSNPdata(y) |> as_tibble()
#names(cousinsDemo)[8:9] = c("ID1", "ID2")

g = getGenotypes(y, ids = ids)

# Dataset with annotation and genotypes
cousinsDem = tibble(FORCE, ID1 = g[1, ], ID2 = g[2, ])

# Add centiMorgan positions
cm = convertPos(chrom = FORCE$CHROM, Mb = FORCE$MB)
cousinsDemo = cousinsDemo |> add_column(CM = cm, .before = "A1")

# Save dataset
usethis::use_data(cousinsDemo, overwrite = TRUE)


# Testing -----------------------------------------------------------------

pars = fitHMM(cousinsDemo)
post = ibdPosteriors(cousinsDemo, k1 = pars$k1, a = pars$a)
segs = findSegments(cousinsDemo, k1 = pars$k1, a = pars$a)
plotIBD(post, segs, ncol = 4)

# With true segments
plotIBD(post, segs, ncol = 4) +
  geom_segment(data = trueSegs, col = "blue", linewidth = 1.5,
               aes(x = startCM, xend = endCM, y = .2, yend = .2))

# k1 estimates
kappaIBD(x, ids)         # pedigree
realisedKappa(sim, ids)  # true realised
ibdEstimate(y, ids)      # thompson estimate

