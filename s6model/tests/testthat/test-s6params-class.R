test_that("A Parameters object is correctly initialized", {
  # Default values
  p <- s6params()

  expect_true(is.s6params(p))
  expect_false(is.s6input(p))

  a <- 0.123
  p.a <- s6params(a = a)
  expect_equal(p.a$a, a)

  Winf <- 162534
  p.Winf <- s6params(Winf = Winf)
  expect_equal(p.Winf$Winf, Winf)

  n <- 0.7
  p.n <- s6params(n = n)
  expect_equal(p.n$n, n)

  expect_warning(s6params(Winf = 1000, Wfs = 1000))
  expect_warning(s6params(Winf = 1000, Wfs = 1010))

  expect_error(s6params(epsilon_a = "a"))
  expect_error(s6params(epsilon_r = "a"))
  expect_error(s6params(a = "a"))
  expect_error(s6params(Winf = "a"))
  expect_error(s6params(Fm = "a"))

  expect_error(s6params(epsilon_a = c(1, 6, 7)))
  expect_error(s6params(epsilon_r = c(1, 6, 7)))
  expect_error(s6params(a = c(1, 6, 7)))
  expect_error(s6params(Winf = c(1, 6, 7)))
  expect_error(s6params(Fm = c(1, 6, 7)))

})


test_that("Initalization works with parameters from local environments", {
  x <- 1.2345
  expect_error(p <- local({
    FmVariable <- x
    s6params(Fm = FmVariable)
  }), NA)
  expect_equal(p$Fm, x)
  expect_equal(p$Fm, 1.2345)
})

test_that("Initialization of Parmaters object using base object", {
  expect_error(s6params(base = "a"))
  expect_error(s6params(base = list(Winf = 1000)))
  expect_error(s6params(base = 1))

  a <- 0.123
  baseobj <- s6params(a = a)

  expect_error(newobj <- s6params(base = baseobj), NA)
  expect_true(all.equal(baseobj, newobj))

  expect_false(difference(baseobj, newobj))

  expect_error(newobj <- s6params(a = 0.654, base = baseobj), NA)
  expect_equal(baseobj$a, 0.123)
  expect_equal(newobj$a, 0.654)
})

test_that("Cloning works as intended", {
  par1 <- s6params(Winf = 1000, Fm = 0.4, Wfs = 100)
  par2 <- par1$.clone()

  expect_false(difference(par1, par2))

  par3 <- par2$.clone()

  expect_true(is.s6params(par1))
  expect_true(is.s6params(par2))
  expect_true(is.s6params(par3))

  expect_true(all.equal(par1, par2))
  expect_true(all.equal(par1, par3))
  expect_true(all.equal(par2, par3))

  expect_false(difference(par1, par3))
  expect_false(difference(par2, par3))
})
