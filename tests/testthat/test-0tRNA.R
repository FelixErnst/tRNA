
context("tRNA")
test_that("tRNA:",{
  expect_error(tRNA:::.checkValueValidity("a", c("b","c")),
               "'a' must be one of the following values: 'b', 'c'")
  
  expect_error(tRNA:::.check_trna_structure_ident("a"),
               "'a' must be one of the following values: ''")
  
})
