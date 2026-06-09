SBIVAR: Spatial BIVARiate association tests across disjoint coordinate
sets<img src='inst/Sbivar.png' align='right' height='20%' width='20%'/>
================

This repo provides code for performing bivariate association tests
between different spatial modalities measured on the same or consecutive
slices, possibly with disjoint coordinate sets. A common coordinate
framework as obtained from alignment is considered given.

<!-- % As introduced in our [preprint](). -->

The package can be installed from GitHub (with suggested packages
through dependencies = TRUE) as follows:

``` r
remotes::install_github("sthawinke/sbivar", build_vignettes = TRUE, dependencies = TRUE)
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

Prepare multithreading

``` r
library(BiocParallel)
nCores <- 2
if (.Platform$OS.type == "unix") {
    # On unix-based systems (linux and macOS), use MulticoreParam
    register(MulticoreParam(nCores))
} else {
    # On windows, use SnowParam
    register(SnowParam(workers = nCores, type = "SOCK"))
}
```

Now analyse this single sample using bivariate Moran’s I.

``` r
moranRes <- sbivar(singleStx, singleMet, singleStxCoords, singleMetCoords,
    method = "Moran's I"
)
```

    ## Starting sbivar analysis of a single image on 2 computing cores

    ## Testing significance of bivariate Moran's I for 25 feature pairs

    ## Calculating bivariate Moran's I statistics ...

    ## Fitting variograms for first modality (5 features) ...

    ## Fitting variograms for second modality (5 features) ...

    ## Calculating variances of bivariate Moran's I statistics ...

Have a look at the results:

``` r
head(moranRes$result)
```

    ##   Modality_X Modality_Y     Ixy_5e.06     Ixy_2e.04      Ixy_0.02 SE.Ixy._5e.06
    ## 1       Pcp4   Dopamine  1.691068e-04  1.313399e-04  6.215548e-05  4.471070e-05
    ## 2    mt.Atp6  Histidine  5.174945e-05  2.628488e-05  6.816206e-06  1.825502e-05
    ## 3       Fth1  Histidine -1.456025e-05 -2.624468e-05 -1.000125e-05  1.823251e-05
    ## 4       Pcp4  Histidine  2.981628e-05  1.236840e-05  3.038781e-06  1.823146e-05
    ## 5       Fth1 Tocopherol -1.413589e-05 -1.422649e-05 -3.450981e-06  1.812031e-05
    ## 6    mt.Atp6    Taurine  1.486419e-06  1.480923e-05  4.817983e-06  1.812031e-05
    ##   SE.Ixy._2e.04 SE.Ixy._0.02         pVal       pAdj
    ## 1  4.147731e-05 2.876564e-05 0.0004216507 0.01054127
    ## 2  7.455806e-06 3.728807e-06 0.0011548134 0.01443517
    ## 3  7.830945e-06 3.390757e-06 0.0019246076 0.01603840
    ## 4  6.869515e-06 3.390383e-06 0.1167989515 0.43101376
    ## 5  7.565471e-06 2.874446e-06 0.1361841917 0.43101376
    ## 6  7.116155e-06 3.146162e-06 0.1462119192 0.43101376

Plot the most significantly spatially associated gene-metabolite pair
Pcp4–Dopamine:

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

    ## Fitting 25 mixed effects models on 2 cores

Extract the results for the desired parameter (the intercept)

``` r
multiGAMLmmsRes <- extractResultsMulti(multiGAMLmms, design)
head(multiGAMLmmsRes$result$Intercept)
```

    ##   Modality_X Modality_Y    Estimate         SE      pVal      pAdj
    ## 1       Pcp4  Histidine -0.13059902 0.07573758 0.1452431 0.9393934
    ## 2       Pcp4 Tocopherol -0.09440682 0.04406096 0.2779911 0.9393934
    ## 3       Gnas    Taurine -0.09714658 0.08080491 0.2830975 0.9393934
    ## 4    mt.Atp6       GABA -0.10968234 0.05900753 0.3142181 0.9393934
    ## 5       Pcp4   Dopamine  0.13036167 0.12926211 0.3594960 0.9393934
    ## 6    mt.Atp6  Histidine  0.11528222 0.12280335 0.3909553 0.9393934

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
