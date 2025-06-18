
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ibdfindr

<!-- badges: start -->

<!-- badges: end -->

The goal of ibdfindr is to detect genomic regions shared *identical by
descent* (IBD) between two individuals, using SNP genotypes. It does so
by fitting a continuous-time hidden Markov model (HMM) to the data.

## Installation

You can install the development version of ibdfindr like so:

``` r
remotes::install_github("magnusdv/ibdfindr")
```

## Example

``` r
library(ibdfindr)
```

As an example we consider the built-in dataset `cousinsDemo`, which
contains SNP genotypes for two related individuals (in fact, a pair of
cousins).

``` r
head(cousinsDemo)
#>   CHROM     MARKER       MB       CM A1 A2  FREQ1 ID1 ID2
#> 1     1  rs9442372 1.083324 0.000000  G  A 0.6110 A/G A/G
#> 2     1  rs4648727 1.844830 0.327833  C  A 0.6372 A/A A/C
#> 3     1 rs10910082 2.487496 1.241475  T  C 0.6260 T/T C/T
#> 4     1  rs6695131 3.084360 1.999883  C  T 0.5719 C/T C/C
#> 5     1  rs3765703 3.675872 4.814675  T  G 0.5942 T/T T/T
#> 6     1  rs7367066 3.887191 5.859379  C  T 0.6739 C/C C/C
```

#### Fitting the HMM

To find IBD segments, we first need to fit a HMM to the data:

``` r
mod = fitHMM(cousinsDemo)
mod
#> $k1
#> [1] 0.1432377
#> 
#> $a
#> [1] 5.025124
```

The `fitHMM()` function returns two parameters: `k1`, the prior
probability of being in an IBD state, and `a`, the switching rate
between IBD/nonIBD states along the chromosome.

#### Finding IBD segments and posteriors

After fitting the model we can find the most likely set of IBD segments
(using the [Viterbi](https://en.wikipedia.org/wiki/Viterbi_algorithm)
algorithm):

``` r
segs = findSegments(cousinsDemo, k1 = mod$k1, a = mod$a)
segs
#>    chrom       start       end   n
#> 1      3  54.8804799 114.69801  68
#> 2      4   0.4048596  41.34093  48
#> 3      4  65.5858206  96.21532  39
#> 4      8 131.9818057 161.53634  39
#> 5      9  27.8427182  60.21339  38
#> 6     12  94.6185515 165.48962  83
#> 7     13   1.3208006  35.75327  44
#> 8     15 100.2667722 117.31752  24
#> 9     16  17.8107159 126.37303 129
#> 10    21   1.7453597  58.62373  68
```

We can also compute the posterior IBD probabilities at each marker locus
(using the
[forward-backward](https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm)
algorithm):

``` r
post = ibdPosteriors(cousinsDemo, k1 = mod$k1, a = mod$a)
head(post)
#>   chrom     marker       cm  freq1  g1  g2 ibs       post
#> 1     1  rs9442372 0.000000 0.6110 1/2 1/2   2 0.07564553
#> 2     1  rs4648727 0.327833 0.6372 2/2 1/2   1 0.07905051
#> 3     1 rs10910082 1.241475 0.6260 1/1 1/2   1 0.08053413
#> 4     1  rs6695131 1.999883 0.5719 1/2 1/1   1 0.08404012
#> 5     1  rs3765703 4.814675 0.5942 1/1 1/1   2 0.09413203
#> 6     1  rs7367066 5.859379 0.6739 1/1 1/1   2 0.08014072
```

The output is a standardised version of the input (in particular, the
genotypes are converted to numeric alleles), with two new columns at the
end: `ibs` (the number of alleles shared IBD) and `post` (the posterior
IBD probability).

#### Plot the results

Finally, we can plot the IBD segments and posteriors. The gray points in
the background show the identity-by-state (IBS) at each marker,
i.e. whether the individuals have 0, 1 or 2 alleles in common. Inferred
IBD regions are shown as red segments at the bottom of each chromosome
panel.

``` r
plotIBD(post, segments = segs, ncol = 4)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />
