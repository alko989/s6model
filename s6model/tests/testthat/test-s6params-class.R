context("Parameters Class")

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


test_that("Initialization of Parmaters object using base object", {
  expect_error(s6params(base = "a"))
  expect_error(s6params(base = list(Winf = 1000)))
  expect_error(s6params(base = 1))
  
  a <- 0.123
  baseobj <- s6params(a = a)
  
  expect_error(newobj <- s6params(base = baseobj), NA)
  expect_true(all.equal(baseobj, newobj))
  
  expect_error(newobj <- s6params(a = 0.654, base = baseobj), NA)
  expect_equal(baseobj$a, 0.123)
  expect_equal(newobj$a, 0.654)
})

