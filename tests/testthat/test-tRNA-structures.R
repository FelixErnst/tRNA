library(tRNA)

context("tRNA structures")
test_that("tRNA structures:",{
  data("gr", package = "tRNA", envir = environment())
  tRNA <- gr
  tRNA <- tRNA[1]
  strList <- gettRNABasePairing(tRNA)
  loopPositions <- .get_loop_positions(strList)
  tRNAStructureTypes <- .get_tRNA_structure_type(loopPositions)
  # discriminator
  discriminator <- tRNA:::.getDiscriminator(tRNA,
                                            strList,
                                            tRNAStructureTypes)
  expect_equal(discriminator, 71)
  # acceptor stem
  acceptor <- tRNA:::.getAcceptorStem(tRNA,
                                      strList,
                                      tRNAStructureTypes,
                                      loopPositions,
                                      discriminator)
  expect_named(acceptor, c("prime5","prime3"))
  expect_named(acceptor[[1]], c("start","end"))
  expect_named(acceptor[[2]], c("start","end"))
  expect_equal(acceptor[[2]], list(start = 64,end = 70))
  # D stem
  dstem <- tRNA:::.getDstem(tRNA,
                            strList,
                            tRNAStructureTypes,
                            loopPositions,
                            acceptor)
  expect_named(dstem, c("prime5","prime3"))
  expect_named(dstem[[1]], c("start","end"))
  expect_named(dstem[[2]], c("start","end"))
  expect_equal(dstem[[2]], list(start = 22,end = 24))
  # T stem
  tstem <- tRNA:::.getTstem(tRNA,
                            strList,
                            tRNAStructureTypes,
                            loopPositions,
                            acceptor)
  expect_named(tstem, c("prime5","prime3"))
  expect_named(tstem[[1]], c("start","end"))
  expect_named(tstem[[2]], c("start","end"))
  expect_equal(tstem[[2]], list(start = 59,end = 63))
  # anticodon stem
  anticodon <- tRNA:::.getAnticodonStem(tRNA,
                                        strList,
                                        tRNAStructureTypes,
                                        loopPositions)
  expect_named(anticodon, c("prime5","prime3"))
  expect_named(anticodon[[1]], c("start","end"))
  expect_named(anticodon[[2]], c("start","end"))
  expect_equal(anticodon[[2]], list(start = 38,end = 42))
  # D loop
  dloop <- tRNA:::.getDloop(tRNA,
                            strList,
                            tRNAStructureTypes,
                            loopPositions,
                            dstem)
  expect_named(dloop, c("start","end"))
  expect_equal(dloop, list(start = 13,end = 21))
  # T loop
  tloop <- tRNA:::.getTloop(tRNA,
                            strList, 
                            tRNAStructureTypes,
                            loopPositions,
                            tstem)
  expect_named(tloop, c("start","end"))
  expect_equal(tloop, list(start = 52,end = 58))
  # anticodon loop
  anticodonloop <- tRNA:::.getAnticodonLoop(tRNA,
                                            strList,
                                            tRNAStructureTypes,
                                            loopPositions,
                                            anticodon)
  expect_named(anticodonloop, c("start","end"))
  expect_equal(anticodonloop, list(start = 31,end = 37))
  # Dprime 5
  dprime5 <- tRNA:::.getDprime5(tRNA,
                                strList,
                                acceptor,
                                dstem)
  expect_named(dprime5, c("start","end"))
  expect_equal(dprime5, list(start = 8,end = 9))
  # Dprime 3
  dprime3 <- tRNA:::.getDprime3(tRNA,
                                strList,
                                dstem,
                                anticodon)
  expect_named(dprime3, c("start","end"))
  expect_equal(dprime3, list(start = 25,end = 25))
  # variable loop
  var <- tRNA:::.getVariableLoop(tRNA,
                                 strList,
                                 tstem,
                                 anticodon)
  expect_named(var, c("start","end"))
  expect_equal(var, list(start = 43,end = 46))
})
