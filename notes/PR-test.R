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


# PR simulations ----------------------------------------------------------

# Number of SNPs
mvec = c(2,4,6,10,20)*1000

prs = sapply(mvec, function(m) {print(m)
  y = x |>
    distributeMarkers(n = m, afreq = c(0.5, 0.5)) |>
    profileSimIBD(sim, verbose = F)

  r = findIBD(y, verbose = F, thomp = F)
  pr = computePR(true7, r$segments)[4:5]
  print(pr)
})
prs
matplot(mvec, t(prs), type = "b")
matplot(
  t(prs), type = "b")
