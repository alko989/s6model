context("fitWL function")

test_that("fitWL chooses the correct columns", {
  a <- 0.00123
  b <- 3.456
  df <- data.frame(bogus = rnorm(51), lebogus = rnorm(51), lengths = 50:100, webogus = rnorm(51), weights =  a * c(50:100) ^ b, weightsbogus = rnorm(51))
  expect_warning(fitab <- fitWL(df, colname.weight = "weights", colname.length = "lengths"), regexp = NA)
  expect_equal(a, fitab$a)
  expect_equal(b, fitab$b)
  expect_equal(nrow(df), fitab$n)
})

test_that("fitWL handles correctly missing values", {
  a <- 0.0234
  b <- 2.987
  df <- data.frame(w = a * c(1,2,3,4,5,0,7,8,9,NA,11,12,13,-1,-1,16) ^ b,
                   l =     c(1,2,3,4,5,6,7,8,9,10,NA,-1,13,-1, 0,NA))
  expect_warning(fitab <- fitWL(df, colname.weight = "w", colname.length = "l", mindata = 8), regexp = NA)
  
  expect_equal(a, fitab$a)
  expect_equal(b, fitab$b)
  expect_equal(9, fitab$n)
  
  expect_warning(defVal <- fitWL(df, colname.weight = "w", colname.length = "l", mindata = 10))
  
  expect_equal(0.01, defVal$a)
  expect_equal(3, defVal$b)
})

test_that("fitWL plots some results", {
  a <- 0.0234
  b <- 2.987
  df <- data.frame(Weight = a * c(1:15) ^ b,
                   Length =     c(1:15))
  fitWL(df, plotFit = FALSE)
  expect_error(recordPlot())
  fitWL(df, plotFit = TRUE)
  expect_error(recordPlot(), NA)
  dev.off()
  expect_warning(fitWL(df, plotFit = TRUE, mindata = 20))
  expect_error(recordPlot(), NA)
  dev.off()
  expect_warning(fitWL(df, plotFit = FALSE, mindata = 20))
  expect_error(recordPlot())
})

test_that("fitWL handles correctly mindata", {
  a <- 0.0123
  b <- 2.345
  smalldf <- data.frame(l = 1:3, w = a * c(1:3) ^ b)
  expect_warning(fit <- fitWL(smalldf, colname.weight = "w", colname.length = "l", mindata = 0))
  expect_warning(fit <- fitWL(smalldf, colname.weight = "w", colname.length = "l", mindata = 1))
  expect_warning(fit <- fitWL(smalldf, colname.weight = "w", colname.length = "l", mindata = 2))
  expect_warning(fit <- fitWL(smalldf, colname.weight = "w", colname.length = "l", mindata = 3))
  expect_warning(fit <- fitWL(smalldf, colname.weight = "w", colname.length = "l", mindata = 4))
  expect_warning(fit <- fitWL(smalldf, colname.weight = "w", colname.length = "l", mindata = 100))
})