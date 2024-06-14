
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ConsensusVelo

<!-- badges: start -->
<!-- badges: end -->

The goal of ConsensusVelo is to estimate RNA velocity with sound
uncertainty quantification in a Bayesian manner. The method consists of
fitting individual genes based on the unspliced and spliced counts,
assuming a time-dependent transcription rate and non-trial initial
conditions. A post-processing step is also provided to combine
information across genes and estimate a gene-shared latent time, which
facilitaes interpretation of cell temporal orders with uncertainty.

A full demo for how to use the R package is provided in
[here](https://github.com/huizizhang949/ConsensusVelo/blob/main/demo/demo.Rmd).

You can install the development version of chdp from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("huizizhang949/ConsensusVelo")
```
