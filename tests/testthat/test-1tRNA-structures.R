
context("tRNA structures")
test_that("tRNA structures:",{
  data("gr", package = "tRNA")
  tRNA <- gr
  tRNA <- tRNA[1]
  #
  strList <- gettRNABasePairing(tRNA)
  expect_equal(colnames(strList[[1]]),c("pos","forward","reverse","character"))
  mcols(tRNA)$tRNA_seq <- as.character(mcols(tRNA)$tRNA_seq)
  expect_equal(strList, gettRNABasePairing(tRNA))
  mcols(tRNA)$tRNA_str <- as.character(mcols(tRNA)$tRNA_str)
  expect_equal(strList, gettRNABasePairing(tRNA))
  tRNA <- gr
  tRNA <- tRNA[1]
  strList <- gettRNABasePairing(tRNA, with.nucleotides = TRUE)
  expect_equal(colnames(strList[[1]]),
               c("pos","forward","reverse","character","base"))
  #
  actual <- gettRNALoopIDs(tRNA)
  expect_s4_class(actual,"LoopIndexList")
  #
  loopPositions <- tRNA:::.get_loop_positions(strList)
  expect_type(loopPositions,"list")
  tRNAStructureTypes <- tRNA:::.get_tRNA_structure_type(loopPositions)
  pos <- lapply(strList,"[[","pos")
  forward <- lapply(strList,"[[","forward")
  reverse <- lapply(strList,"[[","reverse")
  # discriminator
  discriminator <- tRNA:::.getDiscriminator(tRNA,
                                            tRNAStructureTypes)
  expect_equal(discriminator, 71)
  # acceptor stem
  acceptor <- tRNA:::.getAcceptorStem(tRNA,
                                      forward,
                                      reverse,
                                      tRNAStructureTypes,
                                      loopPositions,
                                      discriminator)
  expect_named(acceptor, c("prime5","prime3"))
  expect_named(acceptor[[1]], c("start","end"))
  expect_named(acceptor[[2]], c("start","end"))
  expect_equal(acceptor[[2]], list(start = 64,end = 70))
  # D stem
  dstem <- tRNA:::.getDstem(tRNA,
                            pos,
                            forward,
                            reverse,
                            tRNAStructureTypes,
                            loopPositions,
                            acceptor)
  expect_named(dstem, c("prime5","prime3"))
  expect_named(dstem[[1]], c("start","end"))
  expect_named(dstem[[2]], c("start","end"))
  expect_equal(dstem[[2]], list(start = 22,end = 24))
  # T stem
  tstem <- tRNA:::.getTstem(tRNA,
                            pos,
                            forward,
                            reverse,
                            tRNAStructureTypes,
                            loopPositions,
                            acceptor)
  expect_named(tstem, c("prime5","prime3"))
  expect_named(tstem[[1]], c("start","end"))
  expect_named(tstem[[2]], c("start","end"))
  expect_equal(tstem[[2]], list(start = 59,end = 63))
  # anticodon stem
  anticodon <- tRNA:::.getAnticodonStem(tRNA,
                                        pos,
                                        forward,
                                        reverse,
                                        tRNAStructureTypes,
                                        loopPositions)
  expect_named(anticodon, c("prime5","prime3"))
  expect_named(anticodon[[1]], c("start","end"))
  expect_named(anticodon[[2]], c("start","end"))
  expect_equal(anticodon[[2]], list(start = 38,end = 42))
  # D loop
  dloop <- tRNA:::.getDloop(tRNA,
                            tRNAStructureTypes,
                            loopPositions,
                            dstem)
  expect_named(dloop, c("start","end"))
  expect_equal(dloop, list(start = 13,end = 21))
  # T loop
  tloop <- tRNA:::.getTloop(tRNA,
                            tRNAStructureTypes,
                            loopPositions,
                            tstem)
  expect_named(tloop, c("start","end"))
  expect_equal(tloop, list(start = 52,end = 58))
  # anticodon loop
  anticodonloop <- tRNA:::.getAnticodonLoop(tRNA,
                                            tRNAStructureTypes,
                                            loopPositions,
                                            anticodon)
  expect_named(anticodonloop, c("start","end"))
  expect_equal(anticodonloop, list(start = 31,end = 37))
  # Dprime 5
  dprime5 <- tRNA:::.getDprime5(tRNA,
                                acceptor,
                                dstem)
  expect_named(dprime5, c("start","end"))
  expect_equal(dprime5, list(start = 8,end = 9))
  # Dprime 3
  dprime3 <- tRNA:::.getDprime3(tRNA,
                                dstem,
                                anticodon)
  expect_named(dprime3, c("start","end"))
  expect_equal(dprime3, list(start = 25,end = 25))
  # variable loop
  var <- tRNA:::.getVariableLoop(tRNA,
                                 tstem,
                                 anticodon)
  expect_named(var, c("start","end"))
  expect_equal(var, list(start = 43,end = 46))
  #
  expect_true(hasTloop(tRNA))
  expect_true(hasTStem(tRNA))
  expect_true(hasDStem(tRNA))
  expect_true(hasDloop(tRNA))
  expect_true(hasAcceptorStem(tRNA))
  expect_true(hasAnticodonStem(tRNA))
  expect_true(hasAnticodonLoop(tRNA))
  expect_true(hasVariableLoop(tRNA))
  #
  expect_false(hasTloop(tRNA,length = 10))
  expect_false(hasVariableLoop(tRNA, paired = TRUE))
  expect_false(hasVariableLoop(tRNA, mismatches = TRUE))
  expect_false(hasVariableLoop(tRNA, bulged = TRUE))
  expect_false(hasTStem(tRNA, mismatches = TRUE))
  expect_false(hasTStem(tRNA, bulged = TRUE, mismatches = TRUE))
  expect_false(hasTStem(tRNA, unpaired = TRUE))
})
