library(s6model)

test_that("Empty s6input", {
  inp <- s6input()
  s6.fit(inp)
})
