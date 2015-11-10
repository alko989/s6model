context("fitWL function")

test_that("fitWL chooses the correct columns", {
  a <- 0.00123
  b <- 3.456
  df <- data.frame(bogus = rnorm(51), lebogus = rnorm(51), lengths = 50:100, webogus = rnorm(51), weights =  a * c(50:100) ^ b, weightsbogus = rnorm(51))
  fitab <- fitWL(df, colname.weight = "weights", colname.length = "lengths")
  expect_equal(a, fitab$a)
  expect_equal(b, fitab$b)
  expect_equal(nrow(df), fitab$n)
  
})