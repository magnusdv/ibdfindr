library(pedsuite)
library(ibdsim2)
library(tibble)

# Pedigree
x = cousinPed(1) |> setSNPs(snpData = FORCE)
ids = leaves(x)

# Simulate IBD pattern
sim = ibdsim(x, seed = 123, verbose = FALSE)
realisedKappa(sim, ids)

# True IBD sharing
trueSegs = findPattern(sim, list(carriers = ids))
karyoHaploid(trueSegs)

# Simulate genotypes
y = profileSimIBD(x, sim, ids = ids, seed = 123)

# Quick check of results
# y |> findIBD() |> plotIBD(refSegs = trueSegs)

# Extract and prepare dataset
df = ibdfindr:::getSNPdata(y)
names(df)[8:9] = c("ID1", "ID2")

cousinsDemo = as_tibble(df)

# Save dataset
usethis::use_data(cousinsDemo, overwrite = TRUE)


# Testing -----------------------------------------------------------------

# k1 estimates
kappaIBD(x, ids)              # pedigree
realisedKappa(sim, ids)       # true realised
ibdEstimate(y, ids, verb = F) # thompson estimate
