# test_that("tor_fit return a list of class jagsUI", {
# mod <- tor_fit(MR = c(1,2,1,2,1,2), Ta = c(10,20,10,20,10,20), BMR = 10, TLC= 30, fitting_options = list(ni = 5000, nb = 3000, nc = 1))
# expect_equal(class(mod), "jagsUI")
# })
#
# test_that("tor_fit return an error if MR or Ta are not vectors", {
#   expect_error(tor_fit(MR = data.frame(x = c(1,2,1,2,1,2)), Ta = c(10,20,10,20,10,20), BMR = 10, TLC= 30, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))
# })
#
# test_that("tor_fit return an error if MR or Ta are not the same length", {
#   expect_error(tor_fit(MR = c(1,2,1)), Ta = c(10,20,10,20), BMR = 10, TLC= 30, fitting_options = list(ni = 5000, nb = 3000, nc = 1))
# })
#
# test_that("complete_args completes unclompleted list", {
#
# fun <- function(arg = list(a = 1, b=2, c = 3)){
#   complete_args(fun)
#     out <- arg$c[1]
#   return(out)
# }
# test <- fun(arg = list(c=5))
# expect_equal(test, 5)
# })
