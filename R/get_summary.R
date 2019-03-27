# library(truncnorm)
#
#
# get_summary <- function(mod){
# out <- list()
#
# params <- c("tauy", "inte", "intc", "intr", "betat", "betac", "Tt", "tlc")
#
# mean <- unlist(test$mean[params])
# CI_97.5 <- unlist(test$q97.5[params])
# CI_2.5 <- unlist(test$q2.5[params])
# Rhat <- unlist(test$Rhat[params])
#
# ## frame the output in a df
# x <- as.data.frame(cbind(mean, CI_2.5, CI_97.5, Rhat))
# x$parameter <- rownames(x)
# x <- x[,c(5,1,2,3,4)]
# rownames(x) <- NULL
#
# out$parameter_estimates <- x
# ############### OVERLAP
# prior_tlc <- rtruncnorm(20000,a= test$mean$TLC,mean=0,sd=100)
# prior_Tt <- rtruncnorm(20000,b=test$mean$TLC,mean=0,sd=100)
# prior_betat <- rtruncnorm(20000,b=0,mean=0,sd=100)
#
#
# tlc_overlap <- check_overlap(test, "tlc", prior_tlc)
# Tt_overlap <- check_overlap(test, "Tt", prior_Tt)
# betat_overlap <- check_overlap(test, "betat", prior_betat)
#
# if(Tt_overlap >= 0.3){
#   warning("Parameters Tt and intr not identifiable")
# }
# if(tlc_overlap >= 0.3){
#   warning("Parameters tlc, betac, and intc not identifiable")
# }
# if(betat_overlap >= 0.3){
#   warning("Parameters tlc, betac, and intc not identifiable")
# }
# out$overlap <- list(tlc = tlc_overlap,
#                           Tt = Tt_overlap,
#                           Betat =betat_overlap)
#
# ############## Q10
# out$Q10 <- get_Q102(test)
# test
#
# }
