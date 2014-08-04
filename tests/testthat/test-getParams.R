context("getParams function")

test_that("getParams reads parameters correctly", {

    p1 <- parameters()
    p2 <- p1
    p2@logeta_F <- log(0.01 / p1@scaleeta_F)
    expect_warning(getParams(p2))
    expect_equal(getParams(p2)$eta_F, exp(p1@logeta_F)*p1@scaleeta_F)
    expect_equal(getParams(p2)$p@logeta_F, p1@logeta_F)
    expect_equal(getParams(p2)$p, p1)
    expect_equal(getParams(p1)$p, p1)
})
