context("Plot spot sums")
test_that("Coordinate plotting works", {
    data(Vicari)
    expect_silent(plotSpotSums(
        Vicari$TranscriptOutcomes, Vicari$MetaboliteOutcomes,
        Vicari$TranscriptCoords, Vicari$MetaboliteCoords
    ))
})
