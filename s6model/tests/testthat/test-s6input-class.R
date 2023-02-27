library(s6model)

test_that("An s6input object is correctly initialized", {
  ## Some input data
  set.seed(1)
  data(nscoddat)

  years <- as.numeric(names(nscoddat))
  catch <- data.frame(
    Year = years,
    Catch = runif(length(years), 140, 700))

  inpall <- s6input(wf = nscoddat,
                    catch = catch)

  expect_equal(inpall$wf, nscoddat)
  expect_equal(inpall$years, years)
})

test_that("Subsetting an object works",{
  set.seed(1)
  data(nscoddat)

  years <- as.numeric(names(nscoddat))
  catch <- data.frame(Year = years,
                      Catch = runif(length(nscoddat), 140, 700))


  inpall <- s6input(wf = nscoddat,
                    catch = catch)

  inpten <- inpall[1:10]

  expect_equal(inpten$wf, nscoddat[1:10])
  expect_equal(inpten$wf, nscoddat[1:10])
  expect_equal(inpten$catch, catch[1:10,])
  expect_equal(inpten$years, years[1:10])
  expect_equal(length(inpten), 10)
  expect_equal(seq(inpten), 1:10)
  expect_true(is.s6input(inpall))
  expect_true(is.s6input(inpten))

  expect_error(inpall[1:100])
})

test_that("Printing an object works", {
  set.seed(1)
  data(nscoddat)

  years <- as.numeric(names(nscoddat))
  catch <- data.frame(Year = years,
                      Catch = runif(length(nscoddat), 140, 700))

  inpall <- s6input(wf = nscoddat,
                    catch = catch)
  inpten <- inpall[1:10]

  str <- format(inpall)
  numberoflines <- length(strsplit(str, "\n")[[1]])
  expect_equal(numberoflines, 10)

  str <- format(inpten)
  numberoflines <- length(strsplit(str, "\n")[[1]])
  expect_equal(numberoflines, 8)
  blahbluh <- inpall
  expect_output(print(blahbluh), "blahbluh")
  expect_output(print(blahbluh[1:10]), "blahbluh\\[1:10\\]")
})

test_that("Constructor validation", {
  set.seed(1)
  data(nscoddat)

  years <- as.numeric(names(nscoddat))
  catch <- data.frame(
    Year = years,
    Catch = runif(length(years), 140, 700))

  expect_error(s6input(wf = nscoddat, catch = catch[-1, ]))
})

test_that("Survey data work", {
  set.seed(1)
  data(nscoddat)

  years <- as.numeric(names(nscoddat))
  catch <- data.frame(
    Year = years,
    Catch = runif(length(years), 140, 700))

  inpsur <- s6input(surwf = nscoddat, catch = catch)
  inptensur <- inpsur[1:10]
  expect_equal(inptensur$surwf, nscoddat[1:10])
  expect_equal(inptensur$catch, catch[1:10, ])
  expect_equal(inptensur$years, years[1:10])
  expect_equal(length(inptensur), 10)
  expect_equal(seq(inptensur), 1:10)
  expect_true(is.s6input(inptensur))

  expect_error(s6input(wf = nscoddat, catch = catch[-1,]))
})
