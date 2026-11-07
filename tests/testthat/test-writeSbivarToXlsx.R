if (requireNamespace("openxlsx")) {
    test_that("Writing to spreadsheets works as expected", {
        expect_message(writeSbivarToXlsx(resGAMsSingle, file = "tmpFile", sigLevel = 1))
        expect_true(file.exists("tmpFile.xlsx"))
        expect_error(writeSbivarToXlsx(resGAMsSingle, file = "tmpFile", sigLevel = 1))
        expect_message(writeSbivarToXlsx(resGAMsSingle,
            file = "tmpFile", sigLevel = 1,
            overwrite = TRUE
        ))
        expect_warning(writeSbivarToXlsx(resGAMsSingle,
            file = "tmpFile", sigLevel = 0,
            overwrite = TRUE
        ))
        # Multivariate results
        expect_message(writeSbivarToXlsx(resGAMsMulti,
            file = "tmpFile", sigLevel = 1,
            overwrite = TRUE
        ))
        expect_length(openxlsx::getSheetNames("tmpFile.xlsx"), 3L)
        file.remove("tmpFile.xlsx")
    })
}
