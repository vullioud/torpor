#' plot raw data and model output
#'
#'plot the raw values and the fitted model
#'@name fit_and_plot
#'@param data a fitted model from fit_torpor
#'@param model
#'@return a plot

fit_and_plot <- function(data, model_out = NULL, ...) {

  ## fit a model if necessary
  if(is.null(model_out)) {
    out <- fit_torpor(data = data)
  } else {
    out <- model
  }

  ## define the variables

  beta1<- stats::median(out$sims.list$beta1)
  beta2 <- stats::median(out$sims.list$beta2)
  int1 <- stats::median(out$sims.list$int1)
  int2 <- stats::median(out$sims.list$int2)
  int3 <- stats::median(out$sims.list$int3)
  Tmin <- stats::median(out$sims.list$Tmin)
  TLC <- 	stats::median(out$sims.list$tlc)

  Y <- as.numeric(as.character(data[,2]))
  Ta <- as.numeric(as.character(data[,1]))
  ## remove NAs
  da <- cbind(Y,Ta)[!is.na(Y)&!is.na(Ta)&Ta<(TLC-2),]
  Y <-as.numeric(da[,1])
  Ta <-as.numeric(da[,2])

  Tlimup <- max(Ta,na.rm=T)
  Tlimlo <- min(Ta,na.rm=T)
  MRup <- max(Y,na.rm=T)
  MRlo <- min(Y,na.rm=T)
  ylab <- names(data)[2]


  # get the predictions
  pred <- get_prediction(out, seq(Tlimlo,Tlimup,length=100))
  X <-  pred[pred$Group == "Torp", "Temp"]
  Ymeant <- pred[pred$Group == "Torp", "mean_pred"]
  Y975t <- pred[pred$Group == "Torp", "sup_95"]
  Y025t <- pred[pred$Group == "Torp", "inf_95"]
  Ymeann <- pred[pred$Group == "Norm", "mean_pred"]
  Y975n <- pred[pred$Group == "Norm", "sup_95"]
  Y025n <- pred[pred$Group == "Norm", "inf_95"]


  ## plot
  plot(Y~Ta,type="n",frame=F, xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),ylab=ylab)
  points(Y[out$mean$G<1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))]~Ta[out$mean$G<1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))],col="red",pch=19)
  points(Y[out$mean$G>=1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))]~Ta[out$mean$G>=1.5&(out$Rhat$G<1.1|is.na(out$Rhat$G))],col="blue",pch=19)
  points(Y[out$mean$G>=1.5&out$Rhat$G>=1.1]~Ta[out$mean$G>=1.5&out$Rhat$G>=1.1],col="lightblue",pch=19)
  points(Y[out$mean$G<1.5&out$Rhat$G>=1.1]~Ta[out$mean$G<1.5&out$Rhat$G>=1.1],col="orange",pch=19)
  # points(Y[round(out$mean$G)==1&as.factor(level[state])=="TOR"]~Ta[round(out$mean$G)==1&as.factor(level[state])=="TOR"],pch=3)
  # points(Y[round(out$mean$G)==2&as.factor(level[state])=="EUT"]~Ta[round(out$mean$G)==2&as.factor(level[state])=="EUT"],pch=3)

  if(length(out$mean$G[out$mean$G>1.5&out$Rhat$G<=1.1])>0){
    par(new=T)

    plot(Ymeant~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="blue")
    par(new=T)

    plot(Y975t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="blue")
    par(new=T)

    plot(Y025t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="blue")
    }

  quantile(funnorm(25,int1, beta1),c(0.025,0.5,0.975))

  if(length(out$mean$G[out$mean$G<1.5&out$Rhat$G<=1.1])>0){
    par(new=T)

    plot(Ymeann~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="red")
    par(new=T)

    plot(Y975n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="red")
    par(new=T)

    plot(Y025n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="red")}



  legend("topright",c(
    # paste("BMR:",round(BMR,digits=2)),
    paste("Q10:",round((funtorp(TLC,
                                Tmin,
                                int2,
                                int3,
                                beta1,
                                beta2)/
                          funtorp(Tmin,
                                  Tmin,
                                  int2,
                                  int3,
                                  beta1,
                                  beta2))^(10/(TLC-Tmin)),digits=1)),
    paste("Tmin:",round(Tmin,digits=2)),
    # paste("Tlc:",round(Tlc,digits=2)),
    paste("TBMRt:",round(out$mean$tlc,digits=2)),
    paste("int1:",round(int1,digits=2)),
    paste("int2:",round(int2,digits=2)),
    paste("int3:",round(int3,digits=2)),
    paste("beta1:",round(beta1,digits=2)),
    paste("beta2:",round(beta2,digits=2)),
    # paste("Success:",round(success,digits=2)),
    paste("Rhat:",
          round(max(as.numeric(unlist(out$Rhat)),na.rm=T),digits=2),
          "-",
          round(mean(as.numeric(unlist(out$Rhat)),na.rm=T),digits=2),
          sep="")),
    bty="n",cex=0.8)


}
