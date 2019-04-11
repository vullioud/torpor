#'plot raw data and model output
#'
#'The function provides a plot of the MR values over the respective Ta.
#'Measures are presented in different colors depending of the metabolic
#'stage and regression slopes (continuous lines) as well as 95% CI (segmented lines) are presented.
#'For more flexibility the users are advied to use [fit_torpor()] and [get_prediction()] directly.
#'
#'@name fit_and_plot
#'@param mod a fitted model from [fit_torpor()]
#'@param plot_type A character string specifying the type of plot desired. Either "base" or "ggplot"
#'@param ... arguments to fit the model in [fit_torpor()] if no model is provided
#'@param col_torp color for torpor model fit and points
#'@param col_eut color for euthermy model fit and points
#'@param ylab y label
#'@param xlab x label
#'@param pdf logical if a .pdf copy of the plot should be saved
#'@import ggplot2
#'@export
#'@return a base-R plot or a ggplot object
#'@importFrom grDevices dev.off
#'@examples
#'data(test_data3)
#'fit_and_plot(MR = test_data3[,2], Ta = test_data3[,1], BMR = 1.055, TLC = 29,
#'fitting_options = list(nc =1), plot_type = "ggplot")

fit_and_plot <- function(mod = NULL, plot_type = "ggplot",col_torp = "grey", col_eut = "black", ylab = "MR", xlab = "Ta", pdf = FALSE, ...) {
   # browser()
  ## fit a model if necessary
  if(is.null(mod)) {
    out <- fit_torpor(...)
  } else {
    out <- mod
  }

  TLC <- out$sims.list$TLC[1]
  BMR <- out$sims.list$BMR[1]*out$sims.list$Ym[1]


  Y <- MR <- out$data$Y*out$sims.list$Ym[1]
  Ta <- out$data$Ta

  # plot params
  Tlimup <- max(Ta,na.rm=T)
  Tlimlo <- min(Ta,na.rm=T)
  MRup <- max(Y,na.rm=T)
  MRlo <- min(Y,na.rm=T)


  da <- as.data.frame(cbind(Ta, MR))
  da$G <- as.numeric(out$mean$G)
  da$hat <- out$Rhat$G


  # get the predictions
  pred <- get_prediction(out, seq(Tlimlo,Tlimup,length=100))
  # check overlap
  xxx <- check_overlap(out)


  ###### plot GGplot first
  if(plot_type == "ggplot"){
    G <- lwr_95 <- upr_95 <- NULL
    plot <- ggplot2::ggplot(da, ggplot2::aes(x = Ta, y = MR, col = G > 1.5))+
      ggplot2::geom_point() +
      ggplot2::xlim(c(min(da$Ta), max(da$Ta))) +
      # ylim(c(min(da$MR), max(da$MR))) +
      ggplot2::geom_line(data = pred[pred$group == "Norm", ],
                         ggplot2::aes(x = Ta, y = pred),
                         inherit.aes = F,
                         col = col_eut,
                         linetype = 2) +
      ggplot2::geom_ribbon(data = pred[pred$group == "Norm", ],
                           ggplot2::aes(x = Ta, ymin = lwr_95, ymax = upr_95),
                           inherit.aes = F,
                           alpha = 0.2,
                           fill = col_eut,
                           col = NA)+
      ggplot2::geom_line(data = pred[pred$group == "Torp", ],
                         ggplot2::aes(x = Ta, y = pred),
                         inherit.aes = F,
                         col = col_torp,
                         linetype = 2) +
      ggplot2::geom_ribbon(data = pred[pred$group == "Torp", ],
                           ggplot2::aes(x = Ta, ymin = lwr_95, ymax = upr_95),
                           inherit.aes = F,
                           alpha = 0.2,
                           fill = col_torp,
                           col = NA) +
      ggplot2::scale_color_manual("",
                                    labels = c("Euthermy", "Torpor"),
                                    values =  c(col_eut, col_torp)) +
      scale_y_continuous(name = paste(ylab)) +
      scale_x_continuous(name = paste(xlab))
    if (pdf == TRUE){
      ggplot2::ggsave(plot, filename = "plot.pdf", width = 7, units = "cm")
    }

    return(plot)
  } else {
    X <-  pred[pred$group == "Norm", "Ta"]
    Ymeant <- pred[pred$group == "Torp", "pred"]
    Y975t <- pred[pred$group == "Torp", "upr_95"]
    Y025t <- pred[pred$group == "Torp", "lwr_95"]
    Ymeann <- pred[pred$group == "Norm", "pred"]
    Y975n <- pred[pred$group == "Norm", "upr_95"]
    Y025n <- pred[pred$group == "Norm", "lwr_95"]
    ylab <- paste(ylab)
    xlab <- paste(xlab)

    ## plot
    if(pdf == TRUE){
      pdf("plot.pdf")
    }
    graphics::plot(Y~Ta,type="n",frame=F, xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),ylab=ylab, xlab = xlab)
    graphics::points(Y[out$mean$G<1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))]~Ta[out$mean$G<1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))],col= col_eut,pch=19)
    graphics::points(Y[out$mean$G>=1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))]~Ta[out$mean$G>=1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))],col= col_torp,pch=19)
    graphics::points(Y[out$mean$G>=1.5&out$Rhat$G>=1.1]~Ta[out$mean$G>=1.5&out$Rhat$G>=1.1],col= col_torp,pch=19)
    graphics::points(Y[out$mean$G<1.5&out$Rhat$G>=1.1]~Ta[out$mean$G<1.5&out$Rhat$G>=1.1],col= col_eut ,pch=19)


    if(length(out$mean$G[out$mean$G>1.5&out$Rhat$G<=1.1])>0){
      graphics::par(new=T)

      graphics::plot(Ymeant~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col=col_torp)
      graphics::par(new=T)

      graphics::plot(Y975t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col=col_torp)
      graphics::par(new=T)

      graphics::plot(Y025t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col=col_torp)
    }


    if(length(out$mean$G[out$mean$G<1.5&out$Rhat$G<=1.1])>0){
      graphics::par(new=T)

      graphics::plot(Ymeann~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col=col_eut)
      graphics::par(new=T)

      graphics::plot(Y975n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col=col_eut)
      graphics::par(new=T)

      graphics::plot(Y025n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col=col_eut)}

    graphics::legend("topright",c("Euthermy","Torpor"),pch=19, col=c(col_eut,col_torp),bty="n")
    if(pdf == TRUE){
      dev.off()
    }

  }
}

