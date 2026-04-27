SBIVAR: Spatial BIVARiate association tests across disjoint coordinate
sets<img src='inst/Sbivar.png' align='right' height='20%' width='20%'/>
================

This repo provides code for performing bivariate association tests
between different spatial modalities measured on the same or consecutive
slices, possibly with disjoint coordinate sets. A common coordinate
framework as obtained from alignment is considered given.

<!-- % As introduced in our [preprint](). -->

The package can be installed from GitHub as follows:

``` r
remotes::install_github("sthawinke/sbivar", build_vignettes = TRUE)
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
*MetaboliteOutcomes* for the feature measurements. Have a look at the
aligned coordinates:

``` r
par(mfrow = c(2, 3))
plotCoordsMulti(Vicari$TranscriptCoords, Vicari$MetaboliteCoords, cex = 0.1)
```

![](README_files/figure-gfm/plotvicari-1.png)<!-- -->

``` r
par(mfrow = c(1, 1))
```

## Single-image analysis

The Vicari data consist of six images which could be analysed jointly,
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

    ## Starting sbivar analysis of a single image on 2 computing cores

    ## Testing significance of bivariate Moran's I for 100 feature pairs

    ## Calculating bivariate Moran's I statistics ...

    ## Fitting variograms for first modality (10 features) ...

    ## Fitting variograms for second modality (10 features) ...

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

    ##                        Ixy_5e-06     Ixy_2e-04      Ixy_0.02 SE(Ixy)_5e-06
    ## Fth1__X426.13386   -9.958962e-05 -1.318836e-04 -4.699044e-05  4.884710e-05
    ## Fth1__X573.21671    1.116778e-04  1.238469e-04  6.206438e-05  4.682860e-05
    ## mt.Co1__X426.13386  1.246730e-04  9.255068e-05  1.971629e-05  5.275410e-05
    ## Fth1__X537.21448   -8.479378e-05 -8.609003e-05 -5.844255e-05  4.597603e-05
    ## Fth1__X573.23369    2.839504e-05  4.140745e-05  1.660062e-05  2.663718e-05
    ## mt.Co3__X426.13386  1.090300e-04  7.631569e-05  1.983539e-05  5.517525e-05
    ##                    SE(Ixy)_2e-04 SE(Ixy)_0.02        pVal      pAdj
    ## Fth1__X426.13386    4.492210e-05 3.122381e-05 0.009042104 0.4648595
    ## Fth1__X573.21671    4.302256e-05 3.443858e-05 0.009297191 0.4648595
    ## mt.Co1__X426.13386  4.908292e-05 3.874747e-05 0.042207509 0.9942020
    ## Fth1__X537.21448    4.226008e-05 3.725381e-05 0.062735162 0.9942020
    ## Fth1__X573.23369    1.932606e-05 1.610441e-05 0.082076392 0.9942020
    ## mt.Co3__X426.13386  5.170858e-05 4.130072e-05 0.111022300 0.9942020

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

    ## Starting sbivar analysis (GAMs) of 6 images on 2 computing cores

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

    ## Fitting 100 mixed effects models on 2 cores

Extract the results for the desired parameter (the intercept)

``` r
multiGAMLmmsRes <- extractResultsMulti(multiGAMLmms, design)
head(multiGAMLmmsRes$result$Intercept)
```

    ##                         Estimate          SE       pVal      pAdj
    ## mt.Nd2__X576.20502   0.113286238 0.041491150 0.04126204 0.9014878
    ## Fth1__X555.20713    -0.004059947 0.001738176 0.06673187 0.9014878
    ## mt.Nd4__X573.21671   0.064716265 0.032712328 0.10479963 0.9014878
    ## mt.Nd2__X524.20138  -0.142992590 0.080147167 0.13447853 0.9014878
    ## mt.Nd4__X576.20502   0.130552292 0.078507055 0.15721191 0.9014878
    ## mt.Cytb__X523.19829 -0.167542935 0.108243039 0.18233800 0.9014878

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
