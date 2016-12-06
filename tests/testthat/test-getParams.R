context("getParams function")

test_that("Check if getParams reads parameters correctly", {
    p1 <- parameters()
    p2 <- p1
    p2@logeta_F <- log(0.01)
    expect_warning(gpp2 <- getParams(p2))
    expect_equal(gpp2$eta_F, exp(p1@logeta_F))
    expect_equal(gpp2$p@logeta_F, p1@logeta_F)
    expect_equal(gpp2$p, p1)
    expect_equal(getParams(p1)$p, p1)

    expect_warning(pp <- s6model::getParams(p <- parameters(c("Winf","eta_F"), c(0.5, 0.11), FALSE))$p, regexp = NA)
    expect_equal(pp, p)
})
