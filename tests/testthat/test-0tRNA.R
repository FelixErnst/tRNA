
context("tRNA")
test_that("tRNA:",{
  expect_error(tRNA:::.checkValueValidity("a", c("b","c")),
               "'a' must be one of the following values: 'b', 'c'")
})
