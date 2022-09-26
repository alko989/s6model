context("s6input Class")

test_that("An s6input object is correctly initialized", {
  ## Some input data
  set.seed(1)
  data(nscoddat)
  
  catch <- runif(length(nscoddat), 140, 700)
  years <- as.numeric(names(nscoddat))
  
  inpall <- s6input(wf = nscoddat, isSurvey = FALSE, 
                    years = years,
                    catch = catch)
  
  expect_equal(inpall$wf,       nscoddat)
  expect_equal(inpall$isSurvey, FALSE)
  expect_equal(inpall$years,    years)
})

test_that("Subsetting an object works",{
  set.seed(1)
  data(nscoddat)
  
  catch <- runif(length(nscoddat), 140, 700)
  years <- as.numeric(names(nscoddat))
  
  inpall <- s6input(wf = nscoddat, isSurvey = FALSE, 
                    years = years,
                    catch = catch)
  
  inpten <- inpall[1:10]
  
  expect_equal(inpten$wf, nscoddat[1:10])
  expect_equal(inpten$isSurvey, FALSE)
  expect_equal(inpten$wf, nscoddat[1:10])
  expect_equal(inpten$catch, catch[1:10])
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
  
  catch <- runif(length(nscoddat), 140, 700)
  years <- as.numeric(names(nscoddat))
  
  inpall <- s6input(wf = nscoddat, isSurvey = FALSE, 
                    years = years,
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
})

test_that("Constructor validation", {
  set.seed(1)
  data(nscoddat)
  
  catch <- runif(length(nscoddat), 140, 700)
  years <- as.numeric(names(nscoddat))
  
  expect_error(s6input(wf = nscoddat, isSurvey = FALSE, 
                       years = years[-1],
                       catch = catch))
  expect_error(s6input(wf = nscoddat, isSurvey = FALSE, 
                       years = years,
                       catch = catch[-1]))
})

test_that("Survey data work", {
  set.seed(1)
  data(nscoddat)
  
  catch <- runif(length(nscoddat), 140, 700)
  years <- as.numeric(names(nscoddat))
  
  inpsur <- s6input(surwf = nscoddat, isSurvey = TRUE, 
                    years = years,
                    catch = catch)
  inptensur <- inpsur[1:10]
  expect_equal(inptensur$surwf, nscoddat[1:10])
  expect_equal(inptensur$isSurvey, TRUE)
  expect_equal(inptensur$catch, catch[1:10])
  expect_equal(inptensur$years, years[1:10])
  expect_equal(length(inptensur), 10)
  expect_equal(seq(inptensur), 1:10)
  expect_true(is.s6input(inptensur))
  
  expect_error(s6input(surwf = nscoddat, isSurvey = TRUE, 
                       years = years[-1],
                       catch = catch))
  expect_error(s6input(wf = nscoddat, isSurvey = TRUE, 
                       years = years,
                       catch = catch[-1]))
})