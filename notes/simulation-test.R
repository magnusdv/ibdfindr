library(pedsuite)
library(ibdsim2)

# Pedigree
x = avuncularPed()
ids = leaves(x)

# Simulate IBD pattern
sim = ibdsim(x, ids = ids, seed = 123, verbose = FALSE)
realisedKappa(sim)$perSimulation

# True IBD segments
trueSegs = findPattern(sim, list(carriers = ids))
karyoHaploid(trueSegs)

# OR: Only consider segments > 7 cM
# true7 = findPattern(sim, list(carriers = ids), cutoff = 7, unit = "cm")

# Simulate SNP genotypes
y = x |>
  distributeMarkers(n = 4000, afreq = c(0.5, 0.5)) |>
  profileSimIBD(sim, seed = 123)

# OR: FORCE
# y = x |> setSNPs(snpData = FORCE) |> profileSimIBD(sim, seed = 123)

# Analyse and plot
r = findIBD(y)
plotIBD(r, refSegs = trueSegs)
