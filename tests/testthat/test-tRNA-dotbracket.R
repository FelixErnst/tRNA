library(tRNA)

context("dot bracket annotations")
test_that("dot bracket annotations:",{
  db <- c(">>>>>>>....<<<<<<<")
  actual <- getBasePairing(db)[[1]]
  expect_equal(actual$pos, 1:18)
  expect_equal(actual$forward, actual[order(actual$pos, decreasing = TRUE),]$reverse)
  db <- c(">>>>>>>[[[<<<<<<<]]]{{{{}}}}(((())))")
  actual <- getBasePairing(db)[[1]]
  expect_named(actual,c("pos","forward","reverse","chr"))
  expect_equal(sum(actual$forward),sum(actual$reverse))
  expect_equal(sum(actual$forward),sum(actual$reverse))
  expect_equal(sum(actual[1:7,]$forward),sum(actual[11:17,]$reverse))
  expect_equal(sum(actual[8:10,]$forward),sum(actual[18:20,]$reverse))
  expect_equal(sum(actual[21:24,]$forward),sum(actual[25:28,]$reverse))
  expect_equal(sum(actual[29:32,]$forward),sum(actual[33:36,]$reverse))
  #
  db <- c(">>><<<","<<<>>>")
  actual <- getBasePairing(db)
  expect_equal(actual[[1]],actual[[2]])
  #
  db <- c(">>>>>>><<<<<<<<")
  expect_error(getBasePairing(db),
               "Following structures are invalid:")
  db <- c(")))(((")
  expect_error(getBasePairing(db),
               "Following structures are invalid:")
  #
  db <- c(">>><<<","<<<>>>")
  seq <- c("AGCT")
  expect_error(getBasePairing(db,
                              seq),
               "Number of dot bracket strings and sequences do not match")
  seq <- c("AGCT","AGCT")
  expect_error(getBasePairing(db,
                              seq),
               "Length of dot bracket annotation and sequence do not match for")
  seq <- c("AGCTAA","AGCTAA")
  actual <- getBasePairing(db,
                           seq)
  expect_named(actual[[1]],c("pos","forward","reverse","chr","base"))
  #
  # db <- c(">>>>>>>[[[<<<<<<<]]]{{{{}}}}(((())))")
  db <- c(">>>>>>><<<<<<<{{{{}}}}(((())))")
  actual <- getLoopIDs(db)
  expect_type(actual,"list")
  expect_length(actual[[1]],nchar(db))
  expect_equal(min(actual[[1]]),1)
  expect_equal(max(actual[[1]]),15)
  expect_equal(actual[[1]][length(actual[[1]])],12)
})
