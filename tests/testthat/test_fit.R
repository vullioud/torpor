 test_that("tor_fit return a list of class tor_obj", {
   set.seed(1)
   expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))
   expect_equal(class(mod), c( "tor_obj", "list"))
 })

test_that("tor_fit return an error if M or Ta are not vectors", {
   set.seed(1)
   expect_error(tor_fit(MR = data.frame(x = c(1,2,1,2,1,2)), Ta = c(10,20,10,20,10,20), fitting_options = list(ni = 5000, nb = 3000, nc = 1)))
})

test_that("tor_fit return an error if M or Ta are not the same length", {
   set.seed(1)
   expect_error(tor_fit(MR = c(1,2,1)), Ta = c(10,20,10,20), fitting_options = list(ni = 5000, nb = 3000, nc = 1))
 })

test_that(".complete_args completes unclompleted list", {

 fun <- function(arg = list(a = 1, b=2, c = 3)){
   .complete_args(fun)
    out <- arg$c[1]
   return(out)
 }
   test <- fun(arg = list(c=5))
   expect_equal(test, 5)
})

########### test methods
test_that("method, summary return the same as tor_summarise", {
   set.seed(1)
   expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))
   set.seed(1)
   expect_warning(ref <- tor_summarise(mod))
   set.seed(1)
   expect_warning(job <- summary(mod))
   expect_equal(ref, job)
})

test_that("method, plot return the same as tor_plot", {
   set.seed(1)
   expect_warning(mod <- tor_fit(M = test_data$VO2, Ta = test_data$Ta, fitting_options = list(ni = 5000, nb = 3000, nc = 1)))
   expect_warning(ref <- class(tor_plot(mod, plot_type = "ggplot")))
   expect_warning(job <- class(plot(mod, plot_type = "ggplot")))
   expect_equal(str(ref), str(job))
})

