<img src='inst/Sbivar.png' align='centre' height='15%' width='15%'/>SBIVAR:
Spatial BIVARiate association tests across disjoint coordinate sets
================

This repo provides code for performing bivariate association tests
between different spatial modalities measured on the same or consecutive
slices, possibly with disjoint coordinate sets. A common coordinate
framework (CCF) as obtained from alignment is considered given. As
introduced in our [preprint](https://doi.org/10.1101/2025.05.20.654270).
A simple use-case is shown below, more extensive documentation can be
found in the vignette.

The package can be installed from GitHub as follows:

``` r
library(devtools)
install_github("sthawinke/sbivar")
```

Once installed, you can load the package

``` r
library(sbivar)
```
