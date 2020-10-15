 test_that("tor_summary returns a list of two data.frame", {
    set.seed(123)
   expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))
   expect_warning(test <- tor_summarise(tor_obj = mod))
   length(test)
   expect_equal(length(test), 2)
   expect_equal(class(test[[1]]), "data.frame")
   expect_equal(class(test[[2]]), "data.frame")
 })

 test_that("tor_overlap return a data.frame of 3 rows and 2 columns", {
    set.seed(123)
    expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))
   expect_warning(test <- tor_ppo(tor_obj = mod))
   expect_equal(dim(test), c(5,2))
   expect_equal(class(test)[1], "data.frame")
 })

 test_that("tor_ppo and tor_summarise return an error if the input is not a jagsUI object fitted with tor_fit", {
    set.seed(123)
    data <- data.frame(MR = c(12, 22, 44, 33, 66), Ta = c(10, 20, 30, 40, 50))

   mod <- lm(MR ~ Ta, data = data)
   expect_error(tor_overlap(mod))
   expect_error(tor_summarise(mod))
})
#
test_that("tor_overlap and tor_summarise return a warning if tlc_overlap, Tt_overlap and Betat_overla >= 0.3)", {
   set.seed(123)
  expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))
   expect_warning(tor_ppo(mod))
   expect_warning(tor_summarise(mod))
})

 test_that("get_parameters return a data.frame of 14 rows and 6 col", {
    set.seed(123)
   expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))

   test <- get_parameters(mod)
   expect_equal(nrow(test), 14)
   expect_equal(ncol(test), 6)
})



