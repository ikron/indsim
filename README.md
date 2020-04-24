# indsim
Individual based simulations in R

This is work in progress, building some individual based simulations for R.

Quickstart

Installation
```
devtools::install_github("ikron/indsim") #Note that devtools package is required
```
Load the package
```
library(ggplot2) #Required for plotting functions
library(reshape2) #Required for plotting functions
library(cowplot) #Recommended for plotting function
library(indsim) #Load the package itself
```

Example simulation with infinite alleles model
```
results <- indsim.simulate(N = 1000, generations = 300, sel.intensity = 1, init.f = 0, init.n = 0.5, sigma.e = 1, K = 2000, B = 2, opt.pheno = c(rep(1,100), rep(3,200)), density.reg = "sugar", allele.model = "infinite", mu = 10e-4, n.chr = 3, nloc.chr = c(9,9,9), nqtl = 9, n.delres = 9, delres.ef = 0.1, n.neutral = 9)
```

Packgage manual coming (soon?)