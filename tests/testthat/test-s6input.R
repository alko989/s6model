context("s6input Class")

test_that("A s6input object is correctly initialized", {
    inp <- s6input() ## empty input object
    expect_equal(length(inp), 0)
    
    inp <- simulate(s6params(), binsize = 500, ndataset = 10, years = 1991:2000)
    expect_equal(inp@years, 1991:2000)
    expect_error(s6input(wf = inp@wf, years = 1991:1995))
})