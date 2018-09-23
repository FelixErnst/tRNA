library(tRNA)

context("tRNA structure - non-canonical")
test_that("Human anticodon-, D- and T-loop missing:",{
  data("gr_human", package = "tRNA", envir = environment())
  tRNA <- gr_human
  # check that atRNA with a missing anticodon loop produces the expected result
  # tRNA <- tRNA[1]
  # str <- gettRNAstructureGRanges(tRNA)
  # # discriminator
  # expect_equal(start(str$discriminator), 70)
  # expect_equal(end(str$discriminator), 70)
  # # acceptor stem
  # expect_equal(start(str$acceptorStem$prime5), 1)
  # expect_equal(end(str$acceptorStem$prime5), 7)
  # expect_equal(start(str$acceptorStem$prime3), 63)
  # expect_equal(end(str$acceptorStem$prime3), 69)
  # # Dprime 5
  # expect_equal(start(str$Dprime5), 8)
  # expect_equal(end(str$Dprime5), 9)
  # # D stem
  # expect_equal(start(str$DStem$prime5), 10)
  # expect_equal(end(str$DStem$prime5), 13)
  # expect_equal(start(str$DStem$prime3), 22)
  # expect_equal(end(str$DStem$prime3), 25)
  # # D loop
  # expect_equal(start(str$Dloop), 14)
  # expect_equal(end(str$Dloop), 21)
  # # Dprime 3
  # expect_equal(start(str$Dprime3), 26)
  # expect_equal(end(str$Dprime3), 34)
  # # anticodon stem
  # expect_equal(start(str$anticodonStem$prime5), 35)
  # expect_equal(end(str$anticodonStem$prime5), 38)
  # expect_equal(start(str$anticodonStem$prime3), 39)
  # expect_equal(end(str$anticodonStem$prime3), 42)
  # # anticodon loop
  # expect_equal(start(str$anticodonLoop), 1)
  # expect_equal(end(str$anticodonLoop), 0)
  # # variable loop
  # expect_equal(start(str$variableLoop), 43)
  # expect_equal(end(str$variableLoop), 45)
  # # T stem
  # expect_equal(start(str$TStem$prime5), 46)
  # expect_equal(end(str$TStem$prime5), 50)
  # expect_equal(start(str$TStem$prime3), 58)
  # expect_equal(end(str$TStem$prime3), 62)
  # # T loop
  # expect_equal(start(str$Tloop), 51)
  # expect_equal(end(str$Tloop), 57)
})
