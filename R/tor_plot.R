#'plot raw data and model output
#'
#'The function provides a plot of the MR values over the respective Ta.
#'Measures are presented in different colors depending of the metabolic
#'stage and regression slopes (continuous lines) as well as 95% CI (segmented lines) are presented.
#'For more flexibility the users are advied to use [tor_fit()] and [tor_predict()] directly.
#'
#'@name tor_plot
#'@param mod a fitted model from [tor_fit()]
#'@param plot_type A character string specifying the type of plot desired. Either "base" or "ggplot"
#'@param ... arguments to fit the model in [tor_fit()] if no model is provided
#'@param col_torp color for torpor model fit and points
#'@param col_eut color for euthermy model fit and points
#'@param ylab y label
#'@param xlab x label
#'@param pdf logical if a .pdf copy of the plot should be saved
#'@export
#'@return a base-R plot or a ggplot object
#'@importFrom grDevices dev.off
#'@examples
#'data(test_data3)
#'tor_plot(MR = test_data3[,2], Ta = test_data3[,1], BMR = 1.055, TLC = 29,
#'fitting_options = list(nc =1, ni = 5000, nb = 3000), plot_type = "ggplot")

tor_plot <- function(mod = NULL, plot_type = "ggplot",col_torp = "cornflowerblue", col_eut = "coral3", ylab = "M", xlab = "Ta", pdf = FALSE, ...) {
   # browser()
  ## fit a model if necessary
  if(is.null(mod)) {
    out <- tor_fit(...)
  } else {
    out <- mod
  }

  TLC <- out$sims.list$TLC[1]
  BMR <- out$sims.list$BMR[1]*out$sims.list$Ym[1]


  Y <- MR <- out$data$Y*out$sims.list$Ym[1]
  Ta <- out$data$Ta

  # plot params
  Tlimup <- max(Ta,na.rm=TRUE)
  Tlimlo <- min(Ta,na.rm=TRUE)
  MRup <- max(Y,na.rm=TRUE)
  MRlo <- min(Y,na.rm=TRUE)


  da <- as.data.frame(cbind(Ta, MR))
  da$G <- as.numeric(out$mean$G)
  da$hat <- out$Rhat$G


  # get the predictions
  pred <- tor_predict(out, seq(Tlimlo,Tlimup,length=100))
  # check overlap
  xxx <- tor_overlap(out)


  ###### plot ggplot first
  if(plot_type == "ggplot"){
    G <- lwr_95 <- upr_95 <- NULL
    plot <- ggplot2::ggplot(da, ggplot2::aes(x = Ta, y = MR, col = G > 1.5))+
      ggplot2::geom_point() +
      ggplot2::xlim(c(min(da$Ta), max(da$Ta))) +
      ggplot2::geom_line(data = pred[pred$group == "Euthermy", ],
                         ggplot2::aes(x = Ta, y = pred),
                         inherit.aes = FALSE,
                         col = col_eut,
                         linetype = 2) +
      ggplot2::geom_ribbon(data = pred[pred$group == "Euthermy", ],
                           ggplot2::aes(x = Ta, ymin = lwr_95, ymax = upr_95),
                           inherit.aes = FALSE,
                           alpha = 0.2,
                           fill = col_eut,
                           col = NA)+
      ggplot2::geom_line(data = pred[pred$group == "Torpor", ],
                         ggplot2::aes(x = Ta, y = pred),
                         inherit.aes = FALSE,
                         col = col_torp,
                         linetype = 2) +
      ggplot2::geom_ribbon(data = pred[pred$group == "Torpor", ],
                           ggplot2::aes(x = Ta, ymin = lwr_95, ymax = upr_95),
                           inherit.aes = FALSE,
                           alpha = 0.2,
                           fill = col_torp,
                           col = NA) +
      ggplot2::scale_color_manual("",
                                    labels = c("Euthermy", "Torpor"),
                                    values =  c(col_eut, col_torp)) +
      ggplot2::ylab(paste(ylab)) +
      ggplot2::xlab(paste(xlab))
    if (pdf == TRUE){
      ggplot2::ggsave(plot, filename = "plot.pdf", width = 7, units = "cm")
    }

    return(plot)
  } else {
    Ymeant <- pred[pred$group == "Torpor", "pred"]
    Y975t <- pred[pred$group == "Torpor", "upr_95"]
    Y025t <- pred[pred$group == "Torpor", "lwr_95"]
    Ymeann <- pred[pred$group == "Euthermy", "pred"]
    Y975n <- pred[pred$group == "Euthermy", "upr_95"]
    Y025n <- pred[pred$group == "Euthermy", "lwr_95"]
    ylab <- paste(ylab)
    xlab <- paste(xlab)

    ## plot
    if(pdf == TRUE){
      pdf("plot.pdf")
    }
    graphics::plot(Y~Ta,type="n",frame=FALSE, xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),ylab=ylab, xlab = xlab)
    graphics::points(Y[out$mean$G>=1.5]~Ta[out$mean$G>=1.5],col= col_torp,pch=19)
    graphics::points(Y[out$mean$G<1.5]~Ta[out$mean$G<1.5],col= col_eut ,pch=19)


    if(length(out$mean$G[out$mean$G>1.5])>0){
      X <-  pred[pred$group == "Torpor", "Ta"]
      graphics::par(new=TRUE)

      graphics::plot(Ymeant~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_torp)
      graphics::par(new=TRUE)

      graphics::plot(Y975t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_torp)
      graphics::par(new=TRUE)

      graphics::plot(Y025t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_torp)
    }


    if(length(out$mean$G[out$mean$G<1.5])>0){
      X <-  pred[pred$group == "Euthermy", "Ta"]
      graphics::par(new=TRUE)

      graphics::plot(Ymeann~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_eut)
      graphics::par(new=TRUE)

      graphics::plot(Y975n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_eut)
      graphics::par(new=TRUE)

      graphics::plot(Y025n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_eut)}

    graphics::legend("topright",c("Euthermy","Torpor"),pch=19, col=c(col_eut,col_torp),bty="n")
    if(pdf == TRUE){
      dev.off()
    }

  }
}

