library(pedsuite)
library(ibdsim2)

# Pedigree
x = halfCousinPed(2)
ids = leaves(x)

# Simulate IBD pattern
sim = ibdsim(x, ids = ids, seed = 1, verbose = FALSE)
realisedKappa(sim)

# True IBD segments
trueSegs = findPattern(sim, list(carriers = ids))
karyoHaploid(trueSegs)


# PR simulations ----------------------------------------------------------

# Number of SNPs
nvec = c(2,4,6,10,20)*1000

prs = sapply(nvec, function(m) {
  y = x |>
    distributeMarkers(n = n, afreq = c(0.5, 0.5)) |>
    profileSimIBD(sim, verbose = F)

  r = findIBD(y, verbose = F, thomp = F)
  pr = computePR(true7, r$segments)[4:5] |> cbind(n = n)
  print(pr)
})

prs
matplot(t(prs), type = "b")
