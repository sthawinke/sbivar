#' @rdname sbivar
#' @export
#' @inheritParams sbivarSingle
setMethod("sbivar", "matrix", function(X, Y, Cx, Ey, ...) {
    if(!is.matrix(Y) || !is.matrix(Cx) || !(missing(Ey) || is.matrix(Ey))){
        stop("Since X is a matrix or dataframe, Y, Cx and Ey must be so too!",
             if(is.data.frame(Y) || is.data.frame(Cx) || is.data.frame(Ey)){
                 "\nTry converting data frames with as.matrix()"})
    }
    sbivarSingle(X, Y, Cx, Ey, ...)
})
#' @rdname sbivar
#' @export
#' @inheritParams sbivarMulti
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assay
setMethod("sbivar", "list", function(X, Y, Cx, Ey, assayX = NULL, assayY = NULL, ...) {
    if(!is.list(Y)){
        stop("Since X is a list, Y must be so too!")
    }
    if(all(vapply(X, FUN.VALUE = TRUE, inherits, "SpatialExperiment"))){
        Cx = lapply(X, spatialCoords)
        X = lapply(X, function(y) assay(y, assayX))
    }
    if(all(vapply(Y, FUN.VALUE = TRUE, inherits, "SpatialExperiment"))){
        Ey = lapply(Y, spatialCoords)
        Y = lapply(Y, function(y) assay(y, assayY))
    }
    c(sbivarMulti(X, Y, Cx, Ey, ...), "assayX" = assayX, "assayY" = assayY)
})
#' @param assayX,assayY Assay names to be used in the analysis,
#' see \link[SummarizedExperiment]{assay}
#' @param sample_id_x,sample_id_y If provided, these are used to discriminate between different images
#' included in the same SpatialExperiment objects X and Y. By default, they are assumed to be the same for both X and Y
#' @rdname sbivar
#' @export
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assay
setMethod("sbivar", "SpatialExperiment", function(X, Y, assayX, assayY, sample_id_x,
                                                  sample_id_y = sample_id_x, ...) {
    if(!inherits(Y, "SpatialExperiment")){
        stop("Since X is a SpatialExperiment object, Y must be so too!")
    }
    out = if(missing(sample_id)){
        c(sbivar(t(assay(X, assayX)), t(assay(Y, assayY)), spatialCoords(X), spatialCoords(Y)),
          "assayX" = assayX, "assayY" = assayY, ...)
    } else {
        sbivar(splitSpatialExperiment(X, sample_id_x),
               splitSpatialExperiment(Y, sample_id_y), assayX = assayX, assayY = assayY,...)
    }
    return(out)
})
#' @rdname sbivar
#' @param experiments Names of the experiments in X to be used in the analysis
#' @export
#' @importFrom MultiAssayExperiment MultiAssayExperiment
setMethod("sbivar", "MultiAssayExperiment", function(X, experiments, assayX, assayY, ...) {
    stopifnot(is.character(experiments), length(experiments)==2,
              all(experiments %in% names(X)))
    if(!all(id <- vapply(experiments, FUN.VALUE = TRUE,
                         function(nam) inherits(X[[nam]], "SpatialExperiment")))){
        stop("All components of the MultiAssayExperiment provided must be SpatialExperiment objects.
             Components ", experiments[!id], " are not!")
    }
    sbivar(X[[experiments[1]]], X[[experiments[2]]], assayX = assayX, assayY = assayY, ...)
})
