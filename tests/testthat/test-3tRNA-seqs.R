library(tRNA)

context("tRNA sequences")
test_that("tRNA sequences:",{
  data("gr", package = "tRNA")
  tRNA <- gr[1:10]
  seqs <- gettRNAstructureSeqs(tRNA, joinCompletely = TRUE, joinFeatures = FALSE)
  seqsChars <- gsub("-","",as.character(seqs))
  expect_equal(seqsChars,as.character(tRNA$tRNA_seq))
  expect_named(S4Vectors::metadata(seqs),"tRNA_structures")
  #
  expect_warning(gettRNAstructureSeqs(tRNA, joinCompletely = TRUE,
                                      joinFeatures = TRUE),
                 "Both 'joinCompletely' and 'joinFeatures' are set to TRUE")
  #
  seqs <- gettRNAstructureSeqs(tRNA, joinCompletely = FALSE)
  expect_type(seqs,"list")
  expect_named(seqs,c("acceptorStem","Dprime5","DStem","Dloop","Dprime3",
                      "anticodonStem","anticodonLoop","variableLoop","TStem",
                      "Tloop","discriminator"))
  expect_type(seqs$TStem,"list")
  expect_named(seqs$TStem,c("prime5","prime3"))
  #
  seqs <- gettRNAstructureSeqs(tRNA, joinCompletely = FALSE,
                               joinFeatures = TRUE)
  expect_type(seqs,"list")
  expect_named(seqs,c("acceptorStem","Dprime5","DStem","Dloop","Dprime3",
                      "anticodonStem","anticodonLoop","variableLoop","TStem",
                      "Tloop","discriminator"))
  expect_s4_class(seqs$TStem,"RNAStringSet")
  expect_equal(unique(width(seqs$TStem)), 13L)
  #
  seqs <- gettRNAstructureSeqs(tRNA, joinCompletely = FALSE,
                               joinFeatures = TRUE, padSequences = FALSE)
  expect_type(seqs,"list")
  expect_named(seqs,c("acceptorStem","Dprime5","DStem","Dloop","Dprime3",
                      "anticodonStem","anticodonLoop","variableLoop","TStem",
                      "Tloop","discriminator"))
  expect_s4_class(seqs$TStem,"RNAStringSet")
  expect_named(seqs$TStem)
  expect_equal(unique(width(seqs$TStem)), 10L)
})
