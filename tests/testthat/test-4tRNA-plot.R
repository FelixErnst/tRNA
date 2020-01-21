
context("tRNA plots")
test_that("tRNA plots:",{
  data("gr", package = "tRNA")
  data("gr_eco", package = "tRNA")
  grl <- GRangesList(Sce = gr, Eco = gr_eco)
  plots <- gettRNAFeaturePlots(grl)
  expect_true(all(vapply(plots,is,logical(1),"ggplot")))
})
