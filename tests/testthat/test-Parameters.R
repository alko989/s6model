context("Parameters Class")

test_that("A Parameters object is correctly initialized", {
    p <- parameters() ## Default values

    a <- 0.123
    p.a <- parameters("a", log(a / p@scalea))
    expect_equal(exp(p.a@loga)*p@scalea, a)

    Winf <- 162534
    trWinf <- log(Winf / p@scaleWinf)
    p.trWinf <- parameters("Winf", trWinf, transformed=TRUE)
    p.Winf <- parameters("Winf", Winf, transformed=FALSE)
    expect_equal(exp(p.Winf@logWinf)*p@scaleWinf, Winf)
    expect_equal(exp(p.trWinf@logWinf)*p@scaleWinf, Winf)    
    
    
    p.transformed <- parameters(names = c("a", "A", "Winf", "Wfs"),
                                vals = c(log(0.11 / p@scalea), log(3.91 / p@scaleA),
                                    log(2345 / p@scaleWinf), log(654 / p@scaleWfs)),
                                transformed = TRUE)
    p.nottransformed <- parameters(names = c("a", "A", "Winf", "Wfs"),
                                vals = c(0.11, 3.91,2345, 654),
                                transformed = FALSE)
    expect_true(difference(p.transformed, p.nottransformed))

    expect_equal(parameters(c("Winf", "eta_F"), c(23456, 1234/23456), FALSE),
                 parameters(c("Winf", "Wfs"), c(23456,   1234), FALSE))

   

    expect_warning(parameters(c("Winf", "Wfs"), c(1000, 1000), FALSE))
    expect_warning(parameters(c("Winf", "Wfs"), c(1000, 1001), FALSE))
    expect_warning(etaFWfs <- parameters(c("Winf", "eta_F", "Wfs"), c(1000, 0.1, 800), FALSE))
    expect_equal(etaFWfs@logWfs, log(800/ p@scaleWfs))
    expect_equal(etaFWfs@logeta_F, log(800/1000/ p@scaleeta_F))  
})


test_that("Initialization of Parmaters object using base object", {
    p <- parameters()
    base <- parameters("Winf", 1001, FALSE)
    changed <- parameters("eta_F", 0.4, FALSE, base)
    expect_equal(base@logWinf, changed@logWinf)
    expect_equal(base@loga, changed@loga)
    expect_equal(base@logeta_F, p@logeta_F)
    expect_equal(changed@logeta_F, log(0.4 / p@scaleeta_F))
})

test_that("Winf getter and setter", {
    p <- parameters(c("Winf", "eta_F"), c(1000, parameters()@scaleeta_F), FALSE)
    expect_equal(p@logWfs, log(1000 * p@scaleeta_F / p@scaleWfs))
    expect_equal(Winf(p), 1000)
    expect_equal(p@logeta_F, p@logeta_F)
    Winf(p) <- 1234
    expect_equal(Winf(p), 1234)
})

test_that("Initialization only using Winf sets Wfs/eta_F correctly", {
    p <- parameters("Winf", 3000, FALSE)
    expect_equal(p@logWinf, log(3000/p@scaleWinf))
    expect_equal(p@logWfs, log(3000 * p@scaleeta_F / p@scaleWfs))
})
