setMethod(
    "writeSbivarToXlsx",
    "sbivarResults",
    function(x, file, overwrite, digits, sigLevel) {
        stopifnot(
            is.logical(overwrite), is.character(file), is.numeric(digits),
            is.numeric(sigLevel)
        )
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
        res <- if (x@multi) {
            x@result
        } else {
            list("Baseline" = x@result)
        }
        for (nam in names(res)) {
            mat <- res[[nam]]
            mat <- mat[!is.na(mat[, "pAdj"]), , drop = FALSE]
            # Only significant features
            mat <- mat[mat[, "pAdj"] <= sigLevel, , drop = FALSE]
            if (nrow(mat)) { # Only add sheet when significant findings
                # Rounding
                for (i in setdiff(colnames(mat), c("Modality_X", "Modality_Y"))) {
                    mat[, i] <- signif(mat[, i], digits)
                }
                sheetName <- if (nam == "Intercept") "Baseline" else nam
                addWorksheet(wb, sheetName) # Create sheet and write data to it
                writeData(wb,
                    sheet = sheetName, x = data.frame(mat), colNames = TRUE,
                    rowNames = FALSE
                )
            }
        }
        if (length(wb$worksheets)) {
            saveWorkbook(wb, file = file, overwrite = overwrite)
            message(length(getSheetNames(file)), " tabs successfully written to ", file)
        } else {
            warning(
                "No significant features at significance level ", sigLevel,
                " after multiplicity correction!\nNo file was created."
            )
        }
    }
)
