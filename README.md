SBIVAR: Spatial BIVARiate association tests across disjoint coordinate
sets<img src='inst/Sbivar.png' align='right' height='20%' width='20%'/>
================

This repo provides code for performing bivariate association tests
between different spatial modalities measured on the same or consecutive
slices, possibly with disjoint coordinate sets. A common coordinate
framework (CCF) as obtained from alignment is considered given.

<!-- % As introduced in our [preprint](). -->

Before installing the *sbivar* package, an update of the *smoppix*
package from Github (or BioConductor *devel* branch) of the same author
of at least version 1.1.8 may be needed to gain access to previously
unexported functions. The latest version of the *smoppix* can be
installed as:

``` r
devtools::install_github("sthawinke/smoppix")
```

The package can be installed from GitHub as follows:

``` r
devtools::install_github("sthawinke/sbivar", build_vignettes = TRUE)
```

Once installed, you can load the package

``` r
library(sbivar)
```

As example dataset, we use the one by [Vicari et
al. (2024)](https://doi.org/10.1038/s41587-023-01937-y), which includes
replicated measurements of spatial transcriptomics and metabolomics
coprofiled on mouse brain sections. A subset of the data is available in
the *sbivar* package and can be loaded as:

``` r
data(Vicari)
```

This object consists of four lists of length six: *MetaboliteCoords* and
*TranscriptCoords* for the coordinates, and *TranscriptOutcomes* and
*MetaboliteOutcomes* for the feature measurements. A first sanity check
is to look at the alignment of the coordinates:

``` r
par(mfrow = c(2, 3))
plotCoordsMulti(Vicari$TranscriptCoords, Vicari$MetaboliteCoords, cex = 0.2)
```

![](README_files/figure-gfm/plotvicari-1.png)<!-- -->

``` r
par(mfrow = c(1, 1))
```

## Single-image analysis

The Vicari data consist of six images which are best analysed jointly,
but for didactical purposes we also analyse a single image here, the
sample “V11L12-109_A1”.

``` r
singleSample <- "V11L12-109_A1"
singleStxCoords <- Vicari$TranscriptCoords[[singleSample]]
singleStx <- Vicari$TranscriptOutcomes[[singleSample]]
singleMetCoords <- Vicari$MetaboliteCoords[[singleSample]]
singleMet <- Vicari$MetaboliteOutcomes[[singleSample]]
```

Now analyse this single sample using bivariate Moran’s I.

``` r
library(BiocParallel)
register(MulticoreParam(2))
```

``` r
moranRes <- sbivar(singleStx, singleMet, singleStxCoords, singleMetCoords,
  method = "Moran's I"
)
```

    ## Starting sbivar analysis of a single image on 10 computing cores

    ## Testing significance of bivariate Moran's I for 100 feature pairs

    ## Fitting variograms for first modality (10 features) ...

    ## Fitting variograms for second modality (10 features) ...

    ## Calculating bivariate Moran's I statistics ...

    ## Calculating variances of bivariate Moran's I statistics ...

    ## 10% of tests completed

    ## 20% of tests completed

    ## 30% of tests completed

    ## 40% of tests completed

    ## 50% of tests completed

    ## 60% of tests completed

    ## 70% of tests completed

    ## 80% of tests completed

    ## 90% of tests completed

    ## 100% of tests completed

Have a look at the results:

``` r
head(moranRes$result)
```

    ##                      Ixy_2e-05     Ixy_0.002       Ixy_0.2       pVal     pAdj
    ## Fth1__X573.21671  1.194528e-04  1.022411e-04  1.221550e-05 0.02057933 0.997663
    ## Fth1__X426.13386 -1.185796e-04 -9.845418e-05 -7.282392e-06 0.02677542 0.997663
    ## Fth1__X537.21448 -7.912325e-05 -8.095472e-05 -1.562551e-05 0.11938969 0.997663
    ## Fth1__X573.23369  3.945623e-05  3.128742e-05  3.549317e-06 0.14467783 0.997663
    ## Fth1__X576.20502  7.059640e-05  7.222983e-05  1.586494e-05 0.38430998 0.997663
    ## Fth1__X523.19829 -3.460037e-05 -4.901948e-05 -1.281152e-05 0.43401716 0.997663

Plot the most significantly spatially associated gene-metabolite pair :

``` r
plotTopPair(moranRes, singleStx, singleMet, singleStxCoords, singleMetCoords)
```

![](README_files/figure-gfm/toppairsingle-1.png)<!-- -->

Write the results to a spreadsheet

``` r
writeSbivarToXlsx(moranRes, file = "myfile.xlsx")
```

## Multi-image analysis

Next, we analyse the six images jointly. First construct a variable
identifying the mouse, consisting of the first 10 characters of the
names:

``` r
mouse <- substr(names(Vicari$TranscriptOutcomes), 1, 10)
```

For the multi-image case, we use GAMs as measure of spatial association,
with a negative binomial outcome distribution for the transcriptome
data, and a gamma outcome distribution for the metabolome data. A
log-link is used in both cases.

``` r
multiGAMRes <- sbivar(
  X = Vicari$TranscriptOutcomes, Y = Vicari$MetaboliteOutcomes,
  Cx = Vicari$TranscriptCoords, Ey = Vicari$MetaboliteCoords,
  method = "GAM", families = list("X" = mgcv::nb(), "Y" = Gamma(lin = "log"))
)
```

    ## Starting sbivar analysis of 6 images on 10 computing cores

    ## Image 1 of 6

    ## Image 2 of 6

    ## Image 3 of 6

    ## Image 4 of 6

    ## Image 5 of 6

    ## Image 6 of 6

Next we plug the calculated correlations between the spline surfaces
into a linear model, with random effects for the individual mice:

``` r
design <- data.frame("mouse" = mouse)
multiGAMLmms <- fitLinModels(multiGAMRes, design, Formula = ~ (1 | mouse))
```

    ## Fitting 100 mixed effects models on 10 cores

Extract the results for the desired parameter (the intercept)

``` r
multiGAMLmmsRes <- extractResultsMulti(multiGAMLmms, design)
head(multiGAMLmmsRes$result$Intercept)
```

    ##                       Estimate         SE        pVal      pAdj
    ## Fth1__X426.13386    -0.5173665 0.08692546 0.001913878 0.1913878
    ## mt.Atp6__X426.13386  0.3764458 0.07540160 0.004131023 0.2065512
    ## Gm42418__X555.20345  0.2615645 0.06416528 0.009573444 0.3145332
    ## mt.Co3__X288.1186   -0.7161940 0.19572304 0.014606780 0.3145332
    ## mt.Co1__X288.1186   -0.5606903 0.15622991 0.015726658 0.3145332
    ## Gm42418__X555.20713  0.2342598 0.07195066 0.022549405 0.3678948

No features are significantly associated after multiplicity correction.
For illustration, we plot the feature pair with the smallest p-values
nevertheless:

``` r
plotTopPair(multiGAMLmmsRes,
  Xl = Vicari$TranscriptOutcomes, Yl = Vicari$MetaboliteOutcomes,
  Cxl = Vicari$TranscriptCoords, Eyl = Vicari$MetaboliteCoords, size = 0.3
)
```

![](README_files/figure-gfm/topPairMulti-1.png)<!-- -->

A more extensive description of the sbivar functionality can be found in
the vignette, accessed by calling

``` r
browseVignettes("sbivar")
```
