library(pedsuite)
library(ibdsim2)

# Pedigree
x = halfCousinPed(2)
ids = leaves(x)

# Simulate IBD pattern
sim = ibdsim(x, ids = ids, seed = 1, verbose = FALSE)
realisedKappa(sim)$perSimulation

# True IBD segments
trueSegs = findPattern(sim, list(carriers = ids))
karyoHaploid(trueSegs)

# Simulate genotypes
y = x |> distributeMarkers(n = 10000, afreq = c(0.1, 0.9)) |>
  profileSimIBD(sim, seed = 12)

# Extract SNP data as data frame
w = ibdfindr:::getSNPdata(y)

# Analyse and plot
r = findIBD(w)
plotIBD(r, refSegs = trueSegs)



# Bug
x = halfCousinPed(1)
ids = leaves(x)
sim = ibdsim(x, ids = ids, seed = 123, verbose = FALSE)

y = x |> distributeMarkers(n = 100, afreq = c(0.5, 0.5)) |>
  selectMarkers(chrom=8) |> profileSimIBD(sim, seed = 49)
