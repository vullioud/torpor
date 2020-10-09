 test_that("tor_plot return a object of class ggplot if option plot_type = ggplot", {
  expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))
  expect_warning(x <- tor_plot(mod, plot_type = "ggplot"))
   expect_equal(class(x), c("gg", "ggplot"))
})

