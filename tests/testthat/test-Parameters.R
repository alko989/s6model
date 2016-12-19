context("Parameters Class")

test_that("A Parameters object is correctly initialized", {
    p <- parameters() ## Default values

    a <- 0.123
    p.a <- parameters(list(loga = log(a)))
    expect_equal(exp(p.a@loga), a)

    Winf <- 162534
    trWinf <- log(Winf)
    p.trWinf <- parameters(list(logWinf = trWinf))
    p.Winf <- parameters(list(Winf = Winf))
    expect_equal(exp(p.Winf@logWinf), Winf)
    expect_equal(exp(p.trWinf@logWinf), Winf)    

    p1 <- parameters(c(Winf = Winf, a = a))
    expect_equal(exp(p1@logWfs), Winf * exp(p1@logeta_F))
    
    p.transformed <- parameters(pars =    c(loga = log(0.11), logA = log(3.91), logWinf = log(2345), logWfs = log(654)))
    p.nottransformed <- parameters(pars = c(a    =     0.11,     A =     3.91,     Winf =     2345,     Wfs =     654))
    expect_true(difference(p.transformed, p.nottransformed))

    expect_equal(parameters(c(Winf = 23456, eta_F = 1234/23456)),
                 parameters(c(Winf = 23456, Wfs   = 1234)))
    
    expect_warning(parameters(c(Winf = 1000, Wfs = 1000)))
    expect_warning(parameters(c(Winf = 1000, Wfs = 1001)))
    expect_warning(etaFWfs <- parameters(c(Winf = 1000, eta_F = 0.1, Wfs = 800)))
    expect_equal(etaFWfs@logWfs, log(800))
    expect_equal(etaFWfs@logeta_F, log(800/1000))  
})


test_that("Initialization of Parmaters object using base object", {
    p <- parameters()
    base <- parameters(c(Winf = 1001))
    changed <- parameters(c(eta_F = 0.4), base)
    expect_equal(base@logWinf, changed@logWinf)
    expect_equal(base@loga, changed@loga)
    expect_equal(base@logeta_F, p@logeta_F)
    expect_equal(changed@logeta_F, log(0.4))
})

test_that("Winf getter and setter", {
    p <- parameters(c(Winf = 1000, eta_F = 0.05))
    expect_equal(p@logWfs, log(50))
    expect_equal(getWinf(p), 1000)
    Winf(p) <- 1234
    expect_equal(getWinf(p), 1234)
})

test_that("Initialization only using Winf sets Wfs/eta_F correctly", {
    p <- parameters(c(Winf = 3000))
    expect_equal(p@logWinf, log(3000))
    expect_equal(p@logWfs, log(3000 * 0.05))
})
