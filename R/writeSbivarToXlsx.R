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
#' otherwise throws a warning.
#' @export
#' @importFrom openxlsx createWorkbook writeData addWorksheet saveWorkbook getSheetNames
#' @examples
#' example(sbivar, "sbivar")
#' #The significance level is set to 1 here for illustration,
#' #meaning that all feature pairs will be written to the spreadsheet.
#' # Single result
#' writeSbivarToXlsx(resGAMs, file = "tmpFile", sigLevel = 1)
#' file.remove("tmpFile.xlsx")
#' #Multiple results
#' example(fitLinModels, "sbivar")
#' writeSbivarToXlsx(resMoran, file = "tmpFile", sigLevel = 1)
#' file.remove("tmpFile.xlsx")
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
    wb <- createWorkbook()
    res = switch(results$multiplicity, "multi" = results$result,
               "single" = list("Baseline" = results$result))
    for(nam in names(res)){
        mat = res[[nam]]
        mat <- mat[!is.na(mat[, "pAdj"]), , drop = FALSE]
        # Only significant features
        mat <- mat[mat[, "pAdj"] <= sigLevel, , drop = FALSE]
        if(nrow(mat)){#Only add sheet when significant findings
            # Rounding
            for (i in c("pVal", "pAdj")) {
                mat[, i] <- signif(mat[, i], digits)
            }
            for (i in setdiff(colnames(mat), c("pVal", "pAdj"))) {
                mat[, i] <- round(mat[, i], digits)
            }
            sheetName = if(nam=="Intercept") "Baseline" else nam
            addWorksheet(wb, sheetName) # Create sheet and write data to it
            writeData(wb, sheet = sheetName, x = data.frame(mat), colNames = TRUE,
                      rowNames = TRUE)
        }
    }
    if(length(wb$worksheets)){
        saveWorkbook(wb, file = file, overwrite = overwrite)
        message(length(getSheetNames(file)), " tabs successfully written to ", file)
    } else {
        warning("No significant features at significance level ", sigLevel,
        " after multiplicity correction!\nNo file was created.")
    }
}

