# test_that("tor_predict return a data.frame of 5 columns and n rows ", {
#   mod <- tor_fit(MR = c(1,2,1,2,1,2), Ta = c(10,20,10,20,10,20), BMR = 10, TLC= 30, fitting_options = list(ni = 5000, nb = 3000, nc = 1))
#
#   arg <- c(22, 33)
#   x <- tor_predict(mod, arg)
#   expect_equal(class(x), "data.frame")
#   expect_equal(nrow(x), length(arg))
#   expect_equal(ncol(x), 5)
#
# })
#
#
# test_that("tor_predict return an error if ta is not a vector of numeric", {
#   mod <- tor_fit(MR = c(1,2,1,2,1,2), Ta = c(10,20,10,20,10,20), BMR = 10, TLC= 30, fitting_options = list(ni = 5000, nb = 3000, nc = 1))
#
#   arg <- data.frame(x = 22)
#   arg2 <- "test"
#   expect_error(tor_predict(mod, arg))
#   expect_error(tor_predict(mod, arg2))
# })
#
#
# test_that("tor_classify return a data.frame of length = length(data) and 5 column", {
#   MR <-  c(1,2,1,2,1,2)
#   Ta <- c(10,20,10,20,10,20)
# mod <- tor_fit(MR = MR, Ta = Ta, BMR = 10, TLC= 30, fitting_options = list(ni = 5000, nb = 3000, nc = 1))
#
# x <- tor_classify(mod)
# expect_equal(class(x), "data.frame")
# expect_equal(nrow(x), length(MR))
# expect_equal(ncol(x), 5)
#
# })
#
# test_that("tor_predict_fun correctly backtransform the data", {
#   expect_equal(round(tor_predict_fun(1,1,2,2,2,2,1), 3), 14.778)
# })
#
# test_that("heu_predict_fun correctly backtransform the data", {
#   expect_equal(round(eut_predict_fun(1,1,2,2), 3), 6)
# })
#
