#' Write effect sizes and p-values results to an excel worksheet
#'
#' The results of single- or multi-image analysis are written to an excel spreadsheet
#' with different tabs for every parameter, sorted by increasing p-value.
#'
#' @param results The results of linear model fitting
#' @param file The file to write the results to
#' @param overwrite A boolean, should the file be overwritten if it exists already?
#' @param digits An integer, the number of significant digits to retain for the effect size,
#' raw and adjusted p-values
#' @param sigLevel The significance level threshold to use for the adjusted p-values,
#' only features exceeding the threshold are written to the file. Set this parameter to 1 to write all features
#'
#' @details If no feature exceeds the significance threshold for a certain parameter,
#' an empty tab is created. For each fixed effect, a single tab is written.
#' The "baseline" tabs indicate the overall patterns, the other tabs are named after the fixed effects
#' and indicate departure from this baseline depending on this fixed effect
#' @return Returns invisible with a message when writing operation successful,
#' otherwise throws an error.
#' @export
#' #' @importFrom openxlsx createWorkbook writeData addWorksheet saveWorkbook getSheetNames
writeSbivarToXlsx = function(results, file, overwrite = FALSE, digits = 3,
                             sigLevel = 0.05){
    stopifnot(is.logical(overwrite), is.character(file), is.numeric(digits),
              is.numeric(sigLevel))
    if (!grepl("\\.xlsx", file)) {
        message("Adding .xlsx extension to file")
        file <- paste0(file, ".xlsx")
    }
    if (file.exists(file)) {
        if (overwrite) {
            message("Overwriting existing file")
        } else {
            stop("File ", file, " already exists! Set overwrite = TRUE to overwrite")
        }
    }
    pis <- names(obj)
    fixedEffects <- names(obj[[1]]$results$fixedEffects)
    wb <- createWorkbook()
    for (pi in pis) {
        for (effect in c("Intercept", fixedEffects)) {
            mat <- switch(effect,
                          "Intercept" = obj[[pi]]$results$Intercept,
                          obj[[pi]]$results$fixedEffects[[effect]]
            )
            mat <- mat[!is.na(mat[, "pAdj"]), , drop = FALSE]
            # Only significant features
            subMat <- mat[mat[, "pAdj"] < sigLevel, , drop = FALSE]
            # Rounding
            for (i in c("pVal", "pAdj")) {
                subMat[, i] <- signif(subMat[, i], digits)
            }
            for (i in setdiff(colnames(subMat), c("pVal", "pAdj"))) {
                subMat[, i] <- round(subMat[, i], digits)
            }
            for (smallPI in switch(effect,
                                   "Intercept" = c(TRUE, FALSE),
                                   TRUE
            )) {
                subMat2 <- if (effect == "Intercept") {
                    subMat[match.fun(if (smallPI) "<" else ">")(subMat[, "Estimate"], 0.5), , drop = FALSE]
                } else {
                    subMat
                }
                if(nrow(subMat2)){
                    #Only add sheet when significant findings
                    sheetName <- makeSheetName(pi, effect, smallPI)
                    addWorksheet(wb, sheetName) # Create sheet and write data to it
                    writeData(wb, sheet = sheetName, x = data.frame(subMat2), colNames = TRUE, rowNames = TRUE)
                }
            }
        }
    }
    saveWorkbook(wb, file = file, overwrite = overwrite)
    message(length(getSheetNames(file)), " tabs successfully written to ", file)
}
#' @note Sheet names cannot exceed 31 characters
makeSheetName <- function(pi, effect, smallPI = TRUE) {
    keyword <- if (grepl("Pair", pi)) {
        if (effect != "Intercept") {
            "Colocalization"
        } else if (smallPI) {
            "Colocalized"
        } else {
            "Antilocalized"
        }
    } else if (grepl("nn", pi)) {
        if (smallPI) "Aggregation" else "Regularity"
    } else {
        paste(if (smallPI) "Close to" else "Far from", pi)
    }
    if (grepl("Cell", pi)) {
        keyword <- paste(keyword, "in cell")
    }
    out <- paste0(keyword, "_", switch(effect,
                                       "Intercept" = "baseline",
                                       effect
    ))
    substr(out, 1, min(nchar(out), 31))
}
