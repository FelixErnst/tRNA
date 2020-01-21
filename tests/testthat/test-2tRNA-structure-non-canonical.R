
context("tRNA structure - non-canonical")
test_that("Human anticodon-, D- and T-loop missing:",{
  data("gr_human2", package = "tRNA")
  tRNA <- gr_human2[17]
  # check that atRNA with a missing anticodon loop produces the expected result
  str <- gettRNAstructureGRanges(tRNA)
  # discriminator
  expect_equal(start(str$discriminator), 60)
  expect_equal(end(str$discriminator), 60)
  # acceptor stem
  expect_equal(start(str$acceptorStem$prime5), 1)
  expect_equal(end(str$acceptorStem$prime5), 7)
  expect_equal(start(str$acceptorStem$prime3), 53)
  expect_equal(end(str$acceptorStem$prime3), 59)
  # Dprime 5
  expect_equal(start(str$Dprime5), 8)
  expect_equal(end(str$Dprime5), 9)
  # D stem
  expect_equal(start(str$DStem$prime5), 10)
  expect_equal(end(str$DStem$prime5), 13)
  expect_equal(start(str$DStem$prime3), 23)
  expect_equal(end(str$DStem$prime3), 26)
  # D loop
  expect_equal(start(str$Dloop), 14)
  expect_equal(end(str$Dloop), 22)
  # Dprime 3
  expect_equal(start(str$Dprime3), 27)
  expect_equal(end(str$Dprime3), 27)
  # anticodon stem
  expect_equal(start(str$anticodonStem$prime5), 28)
  expect_equal(end(str$anticodonStem$prime5), 32)
  expect_equal(start(str$anticodonStem$prime3), 40)
  expect_equal(end(str$anticodonStem$prime3), 44)
  # anticodon loop
  expect_equal(start(str$anticodonLoop), 33)
  expect_equal(end(str$anticodonLoop), 39)
  # variable loop
  expect_equal(start(str$variableLoop), 1)
  expect_equal(end(str$variableLoop), 0)
  # T stem
  expect_equal(start(str$TStem$prime5), 1)
  expect_equal(end(str$TStem$prime5), 0)
  expect_equal(start(str$TStem$prime3), 1)
  expect_equal(end(str$TStem$prime3), 0)
  # T loop
  expect_equal(start(str$Tloop), 45)
  expect_equal(end(str$Tloop), 52)
  ##############################################################################
  tRNA <- gr_human2[2]
  str <- gettRNAstructureGRanges(tRNA)
  # D loop
  expect_equal(start(str$Dloop), 8)
  expect_equal(end(str$Dloop), 19)
  # anticodon stem
  expect_equal(start(str$anticodonStem$prime5), 20)
  expect_equal(end(str$anticodonStem$prime5), 23)
  expect_equal(start(str$anticodonStem$prime3), 31)
  expect_equal(end(str$anticodonStem$prime3), 34)
  # Dprime 5
  expect_equal(start(str$Dprime5), 1)
  expect_equal(end(str$Dprime5), 0)
  # Dprime 3
  expect_equal(start(str$Dprime3), 1)
  expect_equal(end(str$Dprime3), 0)
  ##############################################################################
  tRNA <- gr_human2[187]
  str <- gettRNAstructureGRanges(tRNA)
  # D stem
  expect_equal(start(str$DStem$prime3), 21)
  expect_equal(end(str$DStem$prime3), 24)
  # anticodon stem
  expect_equal(start(str$anticodonStem$prime5), 1)
  expect_equal(end(str$anticodonStem$prime5), 0)
  expect_equal(start(str$anticodonStem$prime3), 1)
  expect_equal(end(str$anticodonStem$prime3), 0)
  # anticodon loop
  expect_equal(start(str$anticodonLoop), 25)
  expect_equal(end(str$anticodonLoop), 47)
  # variable loop
  expect_equal(start(str$variableLoop), 1)
  expect_equal(end(str$variableLoop), 0)
  # Dprime 3
  expect_equal(start(str$Dprime3), 1)
  expect_equal(end(str$Dprime3), 0)
  # T stem
  expect_equal(start(str$TStem$prime5), 48)
  expect_equal(end(str$TStem$prime5), 52)
})
