test_that("tor_summary returns a list of two data.frame", {
  mod <- tor_fit(MR = c(12, 22, 44, 33, 66), Ta = c(10, 20, 30, 40, 50), BMR = 20, TLC = 30,
                 fitting_options = list(nc = 2, ni = 5000, nb = 3000))
  test <- tor_summarise(mod = mod)
  length(test)
  expect_equal(length(test), 2)
  expect_equal(class(test[[1]]), "data.frame")
  expect_equal(class(test[[2]]), "data.frame")
})

test_that("tor_overlap return a data.frame of 3 rows and 2 columns", {
  mod <- tor_fit(MR = c(12, 22, 44, 33, 66), Ta = c(10, 20, 30, 40, 50), BMR = 20, TLC = 30,
                 fitting_options = list(nc = 2, ni = 5000, nb = 3000))
  test <- tor_overlap(mod = mod)
  expect_equal(dim(test), c(3,2))
  expect_equal(class(test)[1], "data.frame")
})

test_that("tor_overlap and tor_summarise return an error if the input is not a jagsUI object fitted with tor_fit", {
  data <- data.frame(MR = c(12, 22, 44, 33, 66), Ta = c(10, 20, 30, 40, 50))

  mod <- lm(MR ~ Ta, data = data)
  expect_error(tor_overlap(mod))
  expect_error(tor_summarise(mod))
})

test_that("tor_overlap and tor_summarise return a warning if tlc_overlap, Tt_overlap and Betat_overla >= 0.3)", {

  mod <- tor_fit(MR = c(12, 22, 44, 33, 66), Ta = c(10, 20, 30, 40, 50), BMR = 20, TLC = 30,
                 fitting_options = list(nc = 2, ni = 5000, nb = 3000))
  expect_warning(tor_overlap(mod))
  expect_warning(tor_summarise(mod))
})

