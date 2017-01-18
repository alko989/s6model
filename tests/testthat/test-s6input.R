context("s6input Class")

test_that("A s6input object is correctly initialized", {
    inp <- s6input() ## Default values
    expect_equal(length(inp), 0)
})