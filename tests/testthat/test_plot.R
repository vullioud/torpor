# test_that("tor_plot return a object of class ggplot if option plot_type = ggplot and ", {
#   mod <- tor_fit(MR = c(1,2,1,2,1,2), Ta = c(10,20,10,20,10,20), BMR = 10, TLC= 30, fitting_options = list(ni = 5000, nb = 3000, nc = 1))
#   x <- tor_plot(mod, plot_type = "ggplot")
#   expect_equal(class(x), c("gg", "ggplot"))
# })
#
