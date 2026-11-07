test_that("GAM fitting plotting works", {
    expect_s3_class(plotGAMs(X, Y, Cx, Ey, features = c("X1", "Y2")), "ggplot")
    expect_s3_class(plotGAMs(X, Y, Cx, Ey, features = c("X1", "Y2"), Gamm = TRUE), "ggplot")
    expect_s3_class(plotGAMsTopPair(resGAMsSingle, X, Y, Cx = Cx, Ey = Ey), "ggplot")
    expect_s3_class(plotGAMs(Xl, Yl, Cxl, Eyl, features = c("X1", "Y2")), "ggplot")
    expect_s3_class(plotGAMsTopPair(resGAMsMulti, Xl, Yl, Cx = Cxl, Ey = Eyl), "ggplot")
})
