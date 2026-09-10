#' Plot the fitted splines, and the correlation between them
#'
#' Spawns a three-panel plot with the splines fitted for the two variables, plus a visualization
#' of regions with positive and negative correlations between those splines.
#' The splines are refitted using \link{fitGAM}, so no GAM-objects can be provided.
#' This function takes entire outcome matrices X and Y as argument,
#' to be able to account for offsets in refitting the GAMs.
#' plotGAMsTopResults() plots the feature with the smallest p-value in the 'results' object.
#'
#' @note Both spline surfaces are scaled to the [-1,1] range, the same as
#' the correlation has naturally, for legibility.
#'
#' @inheritParams GAMsSingle
#' @param X,Y Matrices of omics measurements, or lists thereof
#' @param Cx,Ey Corresponding coordinate matrices of dimension two, or lists thereof
#' @param scaleFun The scaling function to be applied before plotting
#' @param addTitle A boolean, should a title be plotted
#' @param features The features to plot
#' @param results Result of a call to \link{sbivar} (single-image) or to
#' \link{fitLinModels} (multi-image)
#' @param ... passed onto \link{fitGAM}
#'
#' @returns A ggplot object
#' @export
#'
#' @examples
#' # Single image
#' example(sbivar, "sbivar")
#' plotGAMs(X, Y, Cx, Ey, features = c("X1", "Y2"))
#' plotGAMsTopResults(resGAMs, X, Y, Cx = Cx, Ey = Ey)
#' # Multi image, arbitrary pair
#' data(Vicari)
#' plotGAMs(Vicari$TranscriptOutcomes, Vicari$MetaboliteOutcomes,
#'     Vicari$TranscriptCoords, Vicari$MetaboliteCoords,
#'     features = c("Pcp4", "Dopamine")
#' )
#' @import ggplot2
#' @importFrom smoppix loadBalanceBplapply
#' @order 1
plotGAMs <- function(X, Y, Cx, Ey, features, scaleFun = "scaleMinusOne",
    families = list("X" = gaussian(), "Y" = gaussian()), addTitle = TRUE, normX = c("none", "rel", "log"),
    normY = c("none", "rel", "log"), n_points_grid = 6e2, Gamm = FALSE, correlation = corGaus(form = ~ x + y, nugget = TRUE, value = c(1, 0.25)), ...) {
    stopifnot(
        is.numeric(n_points_grid), all(vapply(families, FUN.VALUE = TRUE, is, "family")),
        all(vapply(features, FUN.VALUE = TRUE, is.character))
    )
    features <- make.names(features)
    normX <- match.arg(normX)
    normY <- match.arg(normY)
    scaleFun <- get(as.character(scaleFun), mode = "function", getNamespace("sbivar"))
    gamDf <- if (multi <- is.list(X)) {
        foo <- checkInputMulti(X, Y, Cx, Ey)
        gamDfs <- lapply(names(X), function(nam) {
            df <- buildGamDf(
                X[[nam]], Y[[nam]], Cx[[nam]], Ey[[nam]], n_points_grid,
                families, features, scaleFun,
                correlation = correlation,
                normX = normX, normY = normY, Gamm = Gamm
            )$df
            df$image <- nam
            df
        })
        do.call(rbind, gamDfs)
    } else {
        foo <- checkInputSingle(X, Y, Cx, Ey)
        df <- buildGamDf(X, Y, Cx, Ey, n_points_grid, families, features, scaleFun,
            normX = normX, normY = normY, correlation = correlation, Gamm = Gamm
        )
        corEst <- df$corEst
        df$df
    }
    ggplot(gamDf, aes(x, y, fill = Value)) +
        geom_tile() +
        xlab("x-coordinate") +
        ylab("y-coordinate") +
        coord_fixed() +
        facet_grid(if (multi) {
            image ~ feature
        } else {
            ~feature
        }) +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_fill_viridis_c(option = "H", name = "") +
        if (addTitle) {
            labs(title = paste(
                "Spline surfaces, and contributions to correlation estimate",
                if (!multi) round(corEst, 3)
            ))
        }
}
#' @export
#' @rdname plotGAMs
#' @order 2
#' @inheritParams plotTopPair
plotGAMsTopResults <- function(results, X, Y, Cx, Ey, topRank = 1,
    parameter = "Intercept", families = results$families, ...) {
    stopifnot(is.numeric(topRank))
    topFeats <- (
        if (results$multi) {
            results$result[[parameter]]
        } else {
            results$result
        })[topRank, c("Modality_X", "Modality_Y")]
    Cx <- getSpatialCoords(X, Cx)
    X <- getX(X, results$assayX)
    Ey <- getSpatialCoords(Y, Ey)
    Y <- getX(Y, results$assayY)
    plotGAMs(
        X = X, Y = Y, features = topFeats, Cx = Cx, Ey = Ey, families = families,
        multi = results$multi, normX = results$normX, normY = results$normY, Gamm = !results$multi && results$Gamm, correlation = results$correlation, ...
    )
}
#' Make a list of offsets
#'
#' @param X The data matrix
#' @param family The distribution family
#' @returns A list of length two with offsets
makeOffset <- function(X, family) {
    out <- if (family$family != "gaussian") {
        libSizes <- rowSums(X)
        switch(family$link,
            "inverse" = 1 / libSizes,
            "log" = log(libSizes)
        )
    }
    return(out)
}
buildGamDf <- function(X, Y, Cx, Ey, n_points_grid, families, features, scaleFun, normX, normY, pseudoCount = 1e-8, ...) {
    if (families[["X"]]$family != "gaussian") {
        X <- X[idX <- (rowSums(X) > 0), ]
        Cx <- Cx[idX, ]
        if (families[["X"]]$family == "Gamma") {
            X <- X + pseudoCount
        }
    }
    if (families[["Y"]]$family != "gaussian") {
        Y <- Y[idY <- (rowSums(Y) > 0), ]
        Ey <- Ey[idY, ]
        if (families[["Y"]]$family == "Gamma") {
            Y <- Y + pseudoCount
        }
    }
    X <- normMat(X, normX)
    Y <- normMat(Y, normY)
    colnames(Cx) <- colnames(Ey) <- c("x", "y")
    newGrid <- buildNewGrid(Cx = Cx, Ey = Ey, n_points_grid = n_points_grid)
    dfx <- data.frame("value" = X[, features[1]], Cx)
    dfx$Offset <- makeOffset(X, families[["X"]])
    modelx <- fitGAM(
        df = dfx, outcome = "value", family = families[["X"]], ...
    )
    dfy <- data.frame("value" = Y[, features[2]], Ey)
    dfy$Offset <- makeOffset(Y, families[["Y"]])
    modely <- fitGAM(
        df = dfy, outcome = "value", family = families[["Y"]], ...
    )
    predx <- vcovPredGam(modelx, newdata = newGrid)
    predy <- vcovPredGam(modely, newdata = newGrid)
    corContr <- (predx$pred - mean(predx$pred)) * (predy$pred - mean(predy$pred))
    corEst <- sum(corContr) / ((length(predx$pred) - 1) * sd(predx$pred) * sd(predy$pred))
    dat <- rbind(
        data.frame(newGrid, Value = scaleFun(predx$pred), feature = "x"),
        data.frame(newGrid, Value = scaleFun(predy$pred), feature = "y"),
        data.frame(newGrid, Value = corContr / max(abs(corContr), na.rm = TRUE), feature = "cor")
    )
    dat$feature <- factor(dat$feature,
        levels = c("x", "y", "cor"),
        labels = c(features[1], features[2], "Correlation"),
        ordered = TRUE
    )
    return(list("df" = dat, "corEst" = corEst))
}
