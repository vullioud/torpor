#'plot raw data and model output
#'
#'The function provides a plot of the MR values over the respective Ta.
#'Measures are presented in different colors depending of the metabolic
#'stage and regression slopes (continuous lines) as well as 95% CI (segmented lines) are presented.
#'For more flexibility the users are advised to use XXX and XXX
#'@name fit_and_plot
#'@param MR a vector with MR
#'@param Ta, a vector with Ta (same length as MR)
#'@param model_out a fitted model from fit torpor
#'@param ... arguments to fit the model in XXX if not provided
#'@param ggplot1 logical if plot should be done in ggplot
#'@export
#'@return a plot
#'@examples
#'\dontrun{
#'data(test_data2)
#'fit_and_plot(MR = test_data2[,2], Ta = test_data2[,1], BMR = 1.49, TLC = 28.8,
#'fitting_options = list(nc =1), ggplot1 = TRUE)
#'}
fit_and_plot <- function(model_out = NULL, ggplot1 = FALSE, MR, Ta,...) {
  # browser()
  ## fit a model if necessary
  if(is.null(model_out)) {
    out <- fit_torpor(MR = MR, Ta = Ta,...)
  } else {
    out <- model_out
  }
  TLC <- out$sims.list$TLC[1]
  BMR <- out$sims.list$BMR[1]

  ## remove NAs / reshufled the data as done for the fit of the model
  set.seed(666)
   da <- as.data.frame(cbind(MR, Ta)[!is.na(MR) & !is.na(Ta) & Ta < (TLC - 2), ])
   da <- da[sample(nrow(da)), ]

   da$G <- as.numeric(out$mean$G)
   da$hat <- out$Rhat$G


  Y <- da[, "MR"]
  Ta <- da[, "Ta"]

  # plot params
  Tlimup <- max(Ta,na.rm=T)
  Tlimlo <- min(Ta,na.rm=T)
  MRup <- max(Y,na.rm=T)
  MRlo <- min(Y,na.rm=T)
  ylab <- "MR"



  # get the predictions
  pred <- get_prediction(out, seq(Tlimlo,Tlimup,length=100))
  X <-  pred[pred$group == "Norm", "Ta"]
  Ymeant <- pred[pred$group == "Torp", "pred"]
  Y975t <- pred[pred$group == "Torp", "upr_95"]
  Y025t <- pred[pred$group == "Torp", "lwr_95"]
  Ymeann <- pred[pred$group == "Norm", "pred"]
  Y975n <- pred[pred$group == "Norm", "upr_95"]
  Y025n <- pred[pred$group == "Norm", "lwr_95"]

  # if(ggplot == TRUE){
  #   ggplot2::ggplot(da, aes(x = Ta, y = MR, col = G > 1.5))+
  #     ggplot2::geom_point() +
  #     xlim(c(min(da$Ta), max(da$Ta))) +
  #     # ylim(c(min(da$MR), max(da$MR))) +
  #     ggplot2::geom_line(data = pred[pred$group == "Norm", ], ggplot2::aes(x = Ta, y = pred), inherit.aes = F, col = "red", linetype = 2) +
  #     ggplot2::geom_ribbon(data = pred[pred$group == "Norm", ], ggplot2::aes(x = Ta, ymin = lwr_95, ymax = upr_95), inherit.aes = F, alpha = 0.2, fill = "red", col = NA)+
  #     ggplot2::geom_line(data = pred[pred$group == "Torp", ], ggplot2::aes(x = Ta, y = pred), inherit.aes = F, col = "blue", linetype = 2) +
  #     ggplot2::geom_ribbon(data = pred[pred$group == "Torp", ], ggplot2::aes(x = Ta, ymin = lwr_95, ymax = upr_95), inherit.aes = F, alpha = 0.2, fill = "blue", col = NA)
  #
  # } else {

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
