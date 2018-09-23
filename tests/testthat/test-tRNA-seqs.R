library(tRNA)

context("tRNA sequences")
test_that("tRNA sequences:",{
  data("gr", package = "tRNA", envir = environment())
  tRNA <- gr
  seqs <- gettRNAstructureSeqs(gr,
                               joinCompletely = TRUE,
                               joinFeatures = FALSE)
  seqsChars <- gsub("-","",as.character(seqs))
  expect_equal(seqsChars,as.character(tRNA$tRNA_seq))
  expect_named(S4Vectors::metadata(seqs),"tRNA_structures")
})
