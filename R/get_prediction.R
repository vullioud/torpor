#' Get Predictions
#'
#'extract prediction from the output of the model
#'@name get_prediction
#'@param model a fitted model from fit_torpor
#'@param Ta a vector of temperatur for which the prediction should be made
#'@return a data frame with predicted values
#'@examples
#'\dontrun{
#'data(test_data2)
#'test_mod <- fit_torpor(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 98,
#'TLC = 28.88,
#'Model = NULL,
#'fitting_options = list(nc = 1))
#'get_prediction(test_mod, 20)
#'}
#'@export
get_prediction <- function(model, Ta){

  beta1<- (model$sims.list$beta1)
  beta2 <- (model$sims.list$beta2)
  int1 <- (model$sims.list$int1)
  int2 <- (model$sims.list$int2)
  int3 <- (model$sims.list$int3)
  Tmin <- (model$sims.list$Tmin)
  tlc <- (model$sims.list$tlc)

  X <- Ta
  Ymeant<- rep(NA,length(X))
  Y975t<- rep(NA,length(X))
  Y025t<- rep(NA,length(X))
  Ymeann<- rep(NA,length(X))
  Y975n <- rep(NA,length(X))
  Y025n <- rep(NA,length(X))

  for(i in 1:length(Ta)) {
    Ymeant[i] <- stats::median(funtorp(X[i],Tmin, int2, int3, beta1, beta2))
    Y975t[i] <- stats::quantile(funtorp(X[i],Tmin, int2, int3, beta1, beta2),0.975)
    Y025t[i] <- stats::quantile(funtorp(X[i],Tmin, int2, int3, beta1, beta2),0.025)
    Ymeann[i] <- stats::median(funnorm(X[i], int1, beta1))
    Y975n[i] <- stats::quantile(funnorm(X[i], int1, beta1),0.975)
    Y025n[i] <- stats::quantile(funnorm(X[i], int1, beta1),0.025)
  }


#### if model$sims.list$G
out1 <- data.frame(Ta = X,
                  group = rep("Torp", length(X)),
                  pred =  Ymeant,
                  upr_95 = Y975t,
                  lwr_95 = Y025t)

out2 <- data.frame(Ta = X,
                   group = rep("Norm", length(X)),
                   pred =  Ymeann,
                   upr_95 = Y975n,
                   lwr_95 = Y025n)


if(length(model$mean$G[model$mean$G>1.5]) == 0){
  out <- out2
} else if (length(model$mean$G[model$mean$G<1.5]) == 0){
out <- out1
  } else {
out <- rbind(out1, out2)
  }

return(out)

}

#' funtorp
#'
#'fit the model for torpor bats
#'@name funtorp
#'@param x a temperature
#'@param Tmin turning point T
#'@param int2 intercept 1
#'@param int3 intercept 2
#'@param beta1 slope 1
#'@param beta2 slope 2
#'@return a metabolic value

funtorp <- function(x, Tmin, int2, int3, beta1, beta2) {
  out <- ifelse(x<Tmin,int3 +beta1*x,int2*exp(beta2*(x)))
  return(out)
}
#' funnorm
#'
#'fit the model for normotermic bats
#'@name funnorm
#'@param x a temperature
#'@param int1 intercept 1
#'@param beta1 slope 1
#'@return a metabolic value
funnorm<- function(x, int1, beta1) {
  out <- int1 + beta1*x
  return(out)
}


