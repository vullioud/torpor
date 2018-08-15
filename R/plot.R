#'plot raw data and model output
#'
#'plot the raw values and the fitted model
#'@name fit_and_plot
#'@param MR a vector with MR
#'@param Ta, a vector with Ta (same length as MR)
#'@param model_out a fitted model from fit torpor
#'@param ... arguments to fit the model in \code{\link{fit_torpor}} if not provided
#'@export
#'@return a plot
#'@examples
#'data(test_data)
#'fit_and_plot(MR = test_data[,2], Ta = test_data[,1], BMR = 98, TLC = 28.88,
#' fitting_options = list(nc =1) )

fit_and_plot <- function(model_out = NULL, MR, Ta,...) {
# browser()
  ## fit a model if necessary
  if(is.null(model_out)) {
    out <- fit_torpor(MR = MR, Ta = Ta,...)
  } else {
    out <- model_out
  }

  ## define the variables

  beta1<- stats::median(out$sims.list$beta1)
  beta2 <- stats::median(out$sims.list$beta2)
  int1 <- stats::median(out$sims.list$int1)
  int2 <- stats::median(out$sims.list$int2)
  int3 <- stats::median(out$sims.list$int3)
  Tmin <- stats::median(out$sims.list$Tmin)
  tlc <- 	stats::median(out$sims.list$tlc)

  TLC <- out$sims.list$TLC[1]
  BMR <- out$sims.list$BMR[1]

  Y <- as.numeric(as.character(MR))
  Ta <- as.numeric(as.character(Ta))

  ## remove NAs
  da <- cbind(Y,Ta)[!is.na(Y)&!is.na(Ta)&Ta<(TLC-2),] ## TLC estimated -> problem
  Y <-as.numeric(da[,1])
  Ta <-as.numeric(da[,2])

  set.seed(666)
  ord <- order(stats::rnorm(length(Y),0,1))
  Y <- Y[ord]
  Ta <- Ta[ord]

  Tlimup <- max(Ta,na.rm=T)
  Tlimlo <- min(Ta,na.rm=T)
  MRup <- max(Y,na.rm=T)
  MRlo <- min(Y,na.rm=T)
  ylab <- "MR"


  # get the predictions
  pred <- get_prediction(out, seq(Tlimlo,Tlimup,length=100))
  X <-  pred[pred$Group == "Norm", "Temp"]
  Ymeant <- pred[pred$Group == "Torp", "mean_pred"]
  Y975t <- pred[pred$Group == "Torp", "sup_95"]
  Y025t <- pred[pred$Group == "Torp", "inf_95"]
  Ymeann <- pred[pred$Group == "Norm", "mean_pred"]
  Y975n <- pred[pred$Group == "Norm", "sup_95"]
  Y025n <- pred[pred$Group == "Norm", "inf_95"]


  ## plot

  graphics::plot(Y~Ta,type="n",frame=F, xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),ylab=ylab)
  graphics::points(Y[out$mean$G<1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))]~Ta[out$mean$G<1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))],col="red",pch=19)
  graphics::points(Y[out$mean$G>=1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))]~Ta[out$mean$G>=1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))],col="blue",pch=19)
  graphics::points(Y[out$mean$G>=1.5&out$Rhat$G>=1.1]~Ta[out$mean$G>=1.5&out$Rhat$G>=1.1],col="blue",pch=19)
  graphics::points(Y[out$mean$G<1.5&out$Rhat$G>=1.1]~Ta[out$mean$G<1.5&out$Rhat$G>=1.1],col="red",pch=19)


  if(length(out$mean$G[out$mean$G>1.5&out$Rhat$G<=1.1])>0){
    graphics::par(new=T)

    graphics::plot(Ymeant~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="blue")
    graphics::par(new=T)

    graphics::plot(Y975t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="blue")
    graphics::par(new=T)

    graphics::plot(Y025t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="blue")
    }


  if(length(out$mean$G[out$mean$G<1.5&out$Rhat$G<=1.1])>0){
    graphics::par(new=T)

    graphics::plot(Ymeann~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="red")
    graphics::par(new=T)

    graphics::plot(Y975n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="red")
    graphics::par(new=T)

    graphics::plot(Y025n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="red")}

    graphics::legend("topright",c("Normothermy","Torpor"),pch=19, col=c("red","blue"),bty="n")

}
