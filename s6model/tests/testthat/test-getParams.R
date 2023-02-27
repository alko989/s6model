library(testthat)
library(s6model)

test_that("getParams reads parameters correctly", {

    p1 <- s6params()
    p2 <- p1$.clone()
    p2$Wfs <- 123
    gpp2 <- getParams(p2)
    expect_equal(gpp2$Wfs, p2$Wfs)
    expect_equal(gpp2$p, p2)
    expect_true(all.equal(gpp2$p, p2))
    expect_false(all.equal(gpp2$p, p1))
    expect_equal(getParams(p1)$p, p1)

expect_warning(pp <- getParams(p <- s6params(Winf = 0.5, Wfs = 0.11))$p, regexp = NA)
    expect_equal(pp, p)

    p3 <- p1$.clone()
    attr(p3, "s6version") <- "s6model_v0.999.1"
    expect_warning(getParams(p3))
})


test_that("getParams reads parameters correctly", {
   p1 <- s6params()
   expect_error(getParams(p1, calcBRPs = TRUE), NA)
})


test_that("Fmsy is the same from getParams and calcFmsy", {
  p1 <- s6params()
  expect_equal(calcFmsy(p1), getParams(p1, calcBRPs = TRUE)$Fmsy, tolerance = 0.0001)
  p2 <- s6params(epsilon_r = 0.01)
  expect_equal(calcFmsy(p2), getParams(p2, calcBRPs = TRUE)$Fmsy, tolerance = 0.0001)
})
