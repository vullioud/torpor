 test_that("tor_predict return a data.frame of 5 columns and n rows ", {
   set.seed(1)
   expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))

   arg <- c(22, 33)
   set.seed(1)
   expect_warning(x <- tor_predict(mod, arg))
   expect_equal(class(x), "data.frame")
   expect_equal(nrow(x), length(arg))
   expect_equal(ncol(x), 5)

})
#
#
test_that("tor_predict return an error if ta is not a vector of numeric", {
  set.seed(1)
  expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))

  arg <- data.frame(x = 22)
  arg2 <- "test"
  expect_error(tor_predict(mod, arg))
  expect_error(tor_predict(mod, arg2))
 })


test_that("tor_classify return a data.frame of length = length(data) and 5 column", {
   set.seed(1)
   expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))
   x <- tor_classify(mod)
   expect_equal(class(x), "data.frame")
   expect_equal(nrow(x), length(mod$data$Y))
   expect_equal(ncol(x), 4)

})
#
test_that("tor_predict_fun correctly backtransform the data", {
   expect_equal(round(tor_predict_fun(1,1,2,2,2,2,1), 3), 14.778)
})
#
 test_that("heu_predict_fun correctly backtransform the data", {
   expect_equal(round(eut_predict_fun(1,1,2,2), 3), 6)
})
#
