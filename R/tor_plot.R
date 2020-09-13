#'Plot raw data and predicted values
#'
#'[tor_plot()] provides a plot of the metabolic rate (MR) values over the respective ambient temperature (Ta).
#'Raw data and predicted values are presented in different colors depending of the metabolic state.
#`Predicted values are represented by continuous and stripped lines for the estimates’ median and 95CI bounds
#`of the posterior distribution, respectively.
#'For more flexibility the users can use [tor_fit()] and [tor_predict()] directly.
#'
#'@name tor_plot
#'@param tor_obj a fitted model from [tor_fit()]
#'@param plot_type A character string specifying the type of plot desired. Either "base" or "ggplot"
#'@param col_torp color for torpor model prediction and points
#'@param col_eut color for euthermy model prediction and points
#'@param col_n color of NTMZ
#'@param ylab y label
#'@param xlab x label
#'@param pdf logical if a .pdf copy of the plot should be saved
#'@export
#'@return a base-R plot or a ggplot object
#'@importFrom grDevices dev.off

tor_plot <- function(tor_obj = NULL,
                     plot_type = "ggplot",
                     col_torp = "cornflowerblue",
                     col_eut = "coral3",
                     col_n = "green",
                     ylab = "MR",
                     xlab = "Ta",
                     pdf = FALSE) {


  ## please check
  group <- measured_Ta <- measured_MR <- classification <- NULL
  ## retrieve values from the model
  TLC <- tor_obj$out_mtnz_tlc$tlc_estimated
  MTNZ <- tor_obj$out_mtnz_tlc$mtnz_estimated


  Y <- MR <- tor_obj$data$Y
  Ta <- tor_obj$data$Ta

  # plotting limits
  Tlimup <- max(Ta,na.rm=TRUE)
  Tlimlo <- min(Ta,na.rm=TRUE)
  MRup <- max(Y,na.rm=TRUE)
  MRlo <- min(Y,na.rm=TRUE)

  ## get the classification and data
  da <- tor_classify(tor_obj)


  ## get the predictions
  pred <- tor_predict(tor_obj, seq(Tlimlo,Tlimup,length=100)) %>%
    dplyr::rename(classification = group)

  ## check overlap ## warning if > 0.3
  #tor_overlap(out)


  ## ggplot
  if(plot_type == "ggplot"){

    G <- lwr_95 <- upr_95 <- NULL ## check

    plot <- ggplot2::ggplot(da) +
      ggplot2::geom_point(ggplot2::aes(x = measured_Ta, y = measured_MR, col = classification, fill = classification))+
      ggplot2::xlim(c(min(da$measured_Ta), max(da$measured_Ta))) +
      ggplot2::geom_line(data = pred,
                         ggplot2::aes(x = Ta, y = pred, col = classification),
                         inherit.aes = FALSE,
                         linetype = 2) +
      ggplot2::geom_ribbon(data = pred,
                           ggplot2::aes(x = Ta, ymin = lwr_95, ymax = upr_95, fill = classification, group = classification),
                           #inherit.aes = FALSE,
                           alpha = 0.2,
                           col = NA)+
      ggplot2::scale_color_manual("", values =  c(col_eut, col_torp, col_n)) +
      ggplot2::scale_fill_manual("", values =  c(col_eut, col_torp, col_n)) +
      ggplot2::ylab(paste(ylab)) +
      ggplot2::xlab(paste(xlab))+
      ggplot2::theme_light()

    if (pdf == TRUE){
      ggplot2::ggsave(plot, filename = "plot.pdf", width = 7, units = "cm")
    }

    return(plot)

  ## base-R plot
  } else {
    # get values
    Ymeant <- pred[pred$classification == "Torpor", "pred"]
    Y975t <- pred[pred$classification == "Torpor", "upr_95"]
    Y025t <- pred[pred$classification == "Torpor", "lwr_95"]

    Ymeann <- pred[pred$classification == "Euthermy", "pred"]
    Y975n <- pred[pred$classification == "Euthermy", "upr_95"]
    Y025n <- pred[pred$classification == "Euthermy", "lwr_95"]

    YmeanM <- pred[pred$classification == "Ntmz", "pred"]
    Y975M <- pred[pred$classification == "Ntmz", "upr_95"]
    Y025M <- pred[pred$classification == "Ntmz", "lwr_95"]

    ylab <- paste(ylab)
    xlab <- paste(xlab)

   ## plot
    if(pdf == TRUE){
      pdf("plot.pdf")
    }

    graphics::plot(da$measured_MR~da$measured_Ta, type="n",frame=FALSE, xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),ylab= ylab, xlab = xlab)
    graphics::points(da$measured_MR[da$classification == "Torpor"] ~ da$measured_Ta[da$classification == "Torpor"],col = col_torp, pch = 19)
    graphics::points(da$measured_MR[da$classification == "Euthermy"] ~ da$measured_Ta[da$classification == "Euthermy"],col = col_eut , pch = 19)
    graphics::points(da$measured_MR[da$classification == "Ntmz"] ~ da$measured_Ta[da$classification == "Ntmz"],col = col_n , pch = 19)


    if(length(da$classification == "Torpor")> 0){
      X <-  pred[pred$classification == "Torpor", "Ta"]
      graphics::par(new=TRUE)

      graphics::plot(Ymeant~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_torp)
      graphics::par(new=TRUE)

      graphics::plot(Y975t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_torp)
      graphics::par(new=TRUE)

      graphics::plot(Y025t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_torp)
    }


    if(length(da$classification == "Euthermy")>0){
      X <-  pred[pred$classification == "Euthermy", "Ta"]
      graphics::par(new=TRUE)

      graphics::plot(Ymeann~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_eut)
      graphics::par(new=TRUE)

      graphics::plot(Y975n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_eut)
      graphics::par(new=TRUE)

      graphics::plot(Y025n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_eut)

      }

    if(length(da$classification == "Ntmz")>0){
      X <-  pred[pred$classification == "Ntmz", "Ta"]
      graphics::par(new=TRUE)

      graphics::plot(YmeanM~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_n)
      graphics::par(new=TRUE)

      graphics::plot(Y975M~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_n)
      graphics::par(new=TRUE)

      graphics::plot(Y025M~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=FALSE,yaxt="n",ylab="",xlab="",col=col_n)

    }
    graphics::legend("topright",c("Euthermy","Torpor"),pch=19, col=c(col_eut,col_torp),bty="n")
    if(pdf == TRUE){
      dev.off()
    }

  }
}

