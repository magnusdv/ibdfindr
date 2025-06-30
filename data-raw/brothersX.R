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
y = profileSimIBD(x, sim, ids = ids, seed = 1)

# Quick check of results
# y |> findIBD() |> plotIBD(refSegs = trueSegs)

# Extract and prepare dataset
df = ibdfindr:::getSNPdata(y)
names(df)[match(ids, names(df))] = c("ID1", "ID2")

brothersX = as_tibble(df)

# Save dataset
usethis::use_data(brothersX, overwrite = TRUE)

