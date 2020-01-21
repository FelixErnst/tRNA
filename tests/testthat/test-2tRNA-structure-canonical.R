
context("tRNA structure - canonical")
test_that("tRNA structure - canonical:",{
  data("gr", package = "tRNA")
  tRNA <- gr
  expect_true(istRNAGRanges(gr))
  mcols(gr)$tRNA_seq <- NULL
  expect_warning(expect_false(istRNAGRanges(gr)),
                 "Input GRanges object does not meet the requirements")
  # check that a canonical tRNA produces the expected result
  tRNA <- tRNA[1]
  str <- gettRNAstructureGRanges(tRNA)
  # discriminator
  expect_equal(start(str$discriminator), 71)
  expect_equal(end(str$discriminator), 71)
  # acceptor stem
  expect_equal(start(str$acceptorStem$prime5), 1)
  expect_equal(end(str$acceptorStem$prime5), 7)
  expect_equal(start(str$acceptorStem$prime3), 64)
  expect_equal(end(str$acceptorStem$prime3), 70)
  # Dprime 5
  expect_equal(start(str$Dprime5), 8)
  expect_equal(end(str$Dprime5), 9)
  # D stem
  expect_equal(start(str$DStem$prime5), 10)
  expect_equal(end(str$DStem$prime5), 12)
  expect_equal(start(str$DStem$prime3), 22)
  expect_equal(end(str$DStem$prime3), 24)
  # D loop
  expect_equal(start(str$Dloop), 13)
  expect_equal(end(str$Dloop), 21)
  # Dprime 3
  expect_equal(start(str$Dprime3), 25)
  expect_equal(end(str$Dprime3), 25)
  # anticodon stem
  expect_equal(start(str$anticodonStem$prime5), 26)
  expect_equal(end(str$anticodonStem$prime5), 30)
  expect_equal(start(str$anticodonStem$prime3), 38)
  expect_equal(end(str$anticodonStem$prime3), 42)
  # anticodon loop
  expect_equal(start(str$anticodonLoop), 31)
  expect_equal(end(str$anticodonLoop), 37)
  # variable loop
  expect_equal(start(str$variableLoop), 43)
  expect_equal(end(str$variableLoop), 46)
  # T stem
  expect_equal(start(str$TStem$prime5), 47)
  expect_equal(end(str$TStem$prime5), 51)
  expect_equal(start(str$TStem$prime3), 59)
  expect_equal(end(str$TStem$prime3), 63)
  # T loop
  expect_equal(start(str$Tloop), 52)
  expect_equal(end(str$Tloop), 58)
})
