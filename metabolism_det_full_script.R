devtools::load_all()

data(test_data)

  library(rjags)
  library(jagsUI)
  library(MASS)
  library(car)
  library(MCMCvis)
  library(truncnorm)


  ######################################################################
  #Model
  #####################################################################

  #####Enter Data
  #########################


  setwd("/Users/nicolas/Dropbox/Projects/UNIL/Torpor/DATA_TechNote")

  DATA <- read.table("/Users/nicolas/Dropbox/Projects/UNIL/Torpor/DATAExtracted_withoutElephantus.txt", sep = "", header = T)
  liste <- list.files(pattern="*.csv")


  DATA <- test_data
  results <-   matrix(NA,ncol=9,nrow=length(liste))
  colnames(results) <- c("tlc",
                         "bmr",
                         "Tt",
                         "betat",
                         "inte",
                         "betac",
                         "intc",
                         "intr",
                         "MRinhibit")

  level=c("TOR","EUT","TNZ")

  ###########start of the loop
data <- DATA
  data$y <- data$VO2

  for(t in  1:length(liste)){

  data <- read.csv(liste[t],header=T)
  ylab <- names(data)[3]

  names(data)[c(2,3,5)] <- c("Ta","Y","state")

  data <- data[!data$state=="UP",]

  data$state <- gsub("TORREG","TOR",data$state)
  data$state <- gsub("TORTNZ","TNZ",data$state)
  data$state <- gsub("TOR","TOR",data$state)
  data$state <- gsub("HYP","TOR",data$state)
  data$state <- gsub("TOR2","TOR",data$state)

  data <- data[!is.na(data$Y)&!is.na(data$Ta),]
  set.seed(666)
  ord <- order(rnorm(length(data$Y),0,1))
  data <- data[ord,]

  Y <- data$Y
  Ta <- data$Ta
  Ym <- mean(data$Y)

  Tlimup <- max(Ta,na.rm=T)
  Tlimlo <- min(Ta,na.rm=T)
  MRup <- max(Y,na.rm=T)
  MRlo <- min(Y,na.rm=T)
  Title=liste[t]

  ########First step: Find the BMR and get the lowest possible tlc value.
  bmr <- bmr.fit(Y,Ta)[1]
  heterosc<- bmr.fit(Y,Ta)[2]

  ########Second step: Find the tlc value and estimate better BMR
  Ta2 <- Ta[Ta<heterosc]
  Y2 <- Y[Ta<heterosc]

  inits_hetero <- function(){
  	list(
  	  tauy=runif(1,0.05,1),
  	  coef=runif(1,1,1.1),
  	p=runif(2),
  	betat=runif(1,-0.1,0),
  	tlc =runif(1,heterosc,max(Ta)),
  	Tt=runif(1),
  	TMR=runif(1,0,bmr*0.8/Ym),
  	intc1=runif(1,0,bmr*0.8/Ym))}

  params_hetero <- c("tlc","betat","inte")

  win.data <- list(Y = Y2/Ym,
                   NbObservations=length(Y2),
                   Ta=Ta2,
                   heterosc= heterosc,
                   BMR=bmr/Ym,
                   Max=max(Ta))

  # MCMC settings
  ni <- 500000
  nt <- 10
  nb <- 300000
  nc <- 3

  out1 <- jagsUI ::jags(win.data,
                       inits_hetero,
                       params_hetero,
                       paste("hetero1.txt"),
                       n.chains = nc,
                       n.thin = nt,
                       n.iter = ni,
                       n.burnin = nb,
                       parallel=T)

  bmr <- mean(Y[Ta>out1$mean$tlc])
  results[t,2] <- bmr
  results[t,1] <- out1$mean$tlc
  tlc <- out1$mean$tlc

  ########Third step: Assign to torpor, euthermy or TNZ.

  inits_hetero <- function(){
  	list(
  	tauy=runif(1,0.05,1),
  	coef=runif(1,1,1.1),
  	p=runif(3),
  	betat=runif(1,-0.1,0),
  	intc1=runif(1,0,bmr*0.8/Ym),
  	Tt=runif(1),
  	TMR=runif(1,0,bmr*0.8/Ym))}

  params_hetero <- c("tau",
                     "G",
                     "p",
                     "inte",
                     "intc",
                     "intr",
                     "betat",
                     "betac",
                     "Tt",
                     "TMR",
                     "tlc")

  win.data <- list(Y = Y/Ym,
                   NbObservations=length(Y),
                   Ta=Ta,
                   BMR=bmr/Ym,
                   tlc=tlc)

  out2 <- jagsUI ::jags(win.data,
                       inits_hetero,
                       params_hetero,
                       paste("hetero2.txt"),
                       n.chains = nc,
                       n.thin = nt,
                       n.iter = ni,
                       n.burnin = nb,
                       parallel=T)

  X=seq(Tlimlo,tlc,length=100)
  Ymeant=rep(NA,100)
  Y975t=rep(NA,100)
  Y025t=rep(NA,100)
  Ymeann=rep(NA,100)
  Y975n=rep(NA,100)
  Y025n=rep(NA,100)


  #Functions for predicted values
  funtorp <- function(x) (ifelse(x<Tt,intr +betat*x,intc*exp(betac*x)))*Ym
  funnorm<- function(x) (inte +betat*x)*Ym

  betat<-out2$sims.list$betat
  betac <-out2$sims.list$betac
  inte <-out2$sims.list$inte
  intc <-out2$sims.list$intc
  intr <-out2$sims.list$intr
  Tt <-out2$sims.list$Tt
  funtorp(10)

  MRinhibit <- median(funtorp(tlc),na.rm=T)/median(funnorm(tlc),na.rm=T)
  Tm <- (log(bmr)-log(intc))*Tt/(log(funtorp(mean(Tt)))-log(intc))

  for(i in 1:100){
  Ymeant[i] <- median(funtorp(X[i]), na.rm=T)
  Y975t[i] <- quantile(funtorp(X[i]),0.975,na.rm=T)
  Y025t[i] <- quantile(funtorp(X[i]),0.025,na.rm=T)
  Ymeann[i] <- median(funnorm(X[i]),na.rm=T)
  Y975n[i] <- quantile(funnorm(X[i]),0.975,na.rm=T)
  Y025n[i] <- quantile(funnorm(X[i]),0.025,na.rm=T)
  }

  G <- out2$mean$G

  #Define uncertain values (G==0)
  G[sqrt((G-round(G))^2)>0.496] <- 0
  G <- round(G)

  #Any Y value higher than Ymeann should be considered as Euthermic and Y value lower than Ymeant should be considered as torpid
  for(i in 1:length(Y)){
    #Find the X value the closest from Ta[i], from below.
    ind <- seq(1:100)[X==Ta[i]-min((Ta[i]-X)[Ta[i]-X>=0])]
    #Only proceed to this step if there are at least one torpid (G==1) value.
    if(length(G[G==1])!=0){
    G[i] <-
      ifelse(Y[i]>=Ymeann[ind]&G[i]!=3,2,G[i])
    G[i] <-
      ifelse(Y[i]<=Ymeant[ind]&G[i]!=3,1,G[i])}
  }

  data$G <- G

  write.table(data,paste(substr(Title,1,nchar(as.character(Title))-4),".txt",sep=""))


  ########fourth step: Calculate functions without uncertain values

   win.data <- list(Y = Y[G!=0&G!=3]/Ym,
                   NbObservations=length(Y[G!=0&G!=3]),
                   Ta=Ta[G!=0&G!=3],
                   BMR=bmr/Ym,
                   tlc=tlc,
                   G=G[G!=0&G!=3])

  out3 <- jagsUI ::jags(win.data,
                       inits_hetero,
                       params_hetero,
                       paste("hetero3.txt"),
                       n.chains = nc,
                       n.thin = nt,
                       n.iter = ni,
                       n.burnin = nb,
                       parallel=T)

  #Check identifiability of intc, Tt and betat
  PRTintc <- rtruncnorm(20000,a=0,b=bmr/Ym,mean=0,sd=sqrt(1000))
  PRTt <- rtruncnorm(20000,b= tlc,mean=0,sd=sqrt(1000))
  overlapintc <- newfunction(out3, params = 'intc', priors = PRTintc, pdf = FALSE)
  overlapTt <- newfunction(out3, params = 'Tt', priors = PRTt, pdf = FALSE)
  PRbeta <- rtruncnorm(20000,b=0,mean=0,sd=sqrt(100))
  overlapbeta <- newfunction(out3, params = 'betat', priors = PRbeta, pdf = FALSE)

  betat<-out3$sims.list$betat
  betac <-out3$sims.list$betac
  inte <-out3$sims.list$inte
  intc <-out3$sims.list$intc
  intr <-out3$sims.list$intr
  Tt <-out3$sims.list$Tt
  MRinhibit <- median(funtorp(tlc),na.rm=T)/median(funnorm(tlc),na.rm=T)
  Tm <- (log(bmr)-log(intc))*Tt/(log(funtorp(mean(Tt)))-log(intc))

  for(i in 1:100){
    Ymeant[i] <- median(funtorp(X[i]),na.rm=T)
    Y975t[i] <- quantile(funtorp(X[i]),0.975,na.rm=T)
    Y025t[i] <- quantile(funtorp(X[i]),0.025,na.rm=T)
    Ymeann[i] <- median(funnorm(X[i]),na.rm=T)
    Y975n[i] <- quantile(funnorm(X[i]),0.975,na.rm=T)
    Y025n[i] <- quantile(funnorm(X[i]),0.025,na.rm=T)
  }

  Tt <- median(Tt)
  betat<- median(betat)
  inte <- median(inte)
  betac <- median(betac)
  intc <- median(intc)
  intr <- median(intr)
  results[t,3] <- ifelse(overlapTt>30,NA,Tt)
  results[t,4] <- ifelse(overlapbeta>30,NA,betat*Ym)
  results[t,5] <- ifelse(overlapbeta>30,NA,inte*Ym)
  results[t,6] <- ifelse(overlapintc> 30,NA, betac)
  results[t,7] <- ifelse(overlapintc>30,NA,intc*Ym)
  results[t,8] <- ifelse(overlapTt>30,NA,intr*Ym)
  results[t,9] <- MRinhibit

  #########plot
  file <- paste(substr(Title,1,nchar(as.character(Title))-4),".pdf",sep="")
  pdf(file)

  plot(Y~Ta,
       type="n",
       frame=F,
       main=Title,
       xlim=c(Tlimlo, Tlimup),
       ylim=c(MRlo, MRup),
       ylab=ylab)

  points(Y[G==1]~Ta[G==1],
         col="grey",
         pch=19)
  points(Y[G==2]~Ta[G==2],
         col="black",
         pch=19)
  points(Y[G==3]~Ta[G==3])

  points(Y[G==0]~Ta[G==0],pch=19,col="red")



  if(length(G[G==1])>0){
  par(new=T)

  plot(Ymeant~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="grey")
  par(new=T)

  plot(Y975t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="grey")
  par(new=T)

  plot(Y025t~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="grey")}

  if(length(G[G==1])>0& overlapintc>=30& overlapTt<30){
  par(new=T)

  plot(Ymeant[X<Tt]~X[X<Tt],xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="grey")
  par(new=T)

  plot(Y975t[X<Tt]~X[X<Tt],xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="grey")
  par(new=T)

  plot(Y025t[X<Tt]~X[X<Tt],xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="grey")}



  if(length(G[G>1])>0){
  par(new=T)

  plot(Ymeann~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="black")
  par(new=T)

  plot(Y975n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="black")
  par(new=T)

  plot(Y025n~X,xlim=c(Tlimlo, Tlimup),ylim=c(MRlo, MRup),type="l",lty=2,xaxt="n",frame=F,yaxt="n",ylab="",xlab="",col="black")

  segments(tlc, bmr,Tlimup,bmr)
  }


  if(length(G[G==1])>0){

  legend("topright",c(
  paste("Tauc:",round(out3$mean$tau[1],digits=3)),
  paste("Taur:",round(out3$mean$tau[2],digits=3)),
  paste("BMR:",round(bmr,digits=3)),
  paste("TMR:",round(out3$mean$TMR*Ym,digits=3)),
  paste("Tlc:",round(tlc,digits=3)),
  paste("Tt:",ifelse(overlapTt> 30,"NI",round(Tt,digits=3))),
  paste("inte:",round(inte*Ym,digits=3)),
  paste("intc:",ifelse(overlapintc> 30,"NI",round(intc*Ym,digits=3))),
  paste("intr:",ifelse(overlapTt> 30,"NI",round(intr*Ym,digits=3))),
  paste("betat:",ifelse(overlapbeta> 30,"NI",round(betat*Ym,digits=3))),
  paste("betac:",ifelse(overlapintc> 30,"NI",round(betac,digits=3))),
  paste("MRinhibit:",round(MRinhibit,3)),
  paste("Rhat:",
  ifelse(
  all(out3$Rhat$tau>=1.1)|all(unlist(out3$Rhat)[(length(unlist(out3$Rhat))-8):(length(unlist(out3$Rhat))-2)]>=1.1),"not converged","OK"),
        sep=" ")),
  bty="n",cex=0.8)
    }else{
  legend("topright",c(
    paste("Tauc:",round(out3$mean$tau[1],digits=3)),
    paste("Taur:",round(out3$mean$tau[2],digits=3)),
  paste("BMR:",round(bmr,digits=3)),
  paste("inte:",round(inte*Ym,digits=3)),
  paste("betat:",ifelse(overlapbeta> 30,"NI",round(betat*Ym,digits=3))),
  paste("Rhat:",
  ifelse(
  all(out3$Rhat$tau>=1.1)|all(unlist(out3$Rhat)[(length(unlist(out3$Rhat))-8):(length(unlist(out3$Rhat))-1)]>=1.1),"not converged","OK"),
        sep=" ")),
  bty="n",cex=0.8)}

  abline(h=DATA[t,4],col="red")
  abline(v=DATA[t,2],col="red")
  abline(v=tlc,col="grey")

  state <- data$state
  state <-gsub("TOR","1",state)
  state <-gsub("EUT","2",state)
  state <-gsub("TNZ","3",state)

  points(Y[G!=state&G!=0]~Ta[G!=state&G!=0],pch=3)


  dev.off( )
  }

  write.table(results,"/Users/nicolas/Dropbox/Projects/UNIL/Torpor/TLCBMR.txt")
  data.results <- read.table("/Users/nicolas/Dropbox/Projects/UNIL/Torpor/TLCBMR.txt")


  plot(1,frame=F,xaxt="n",yaxt="n",xlab=expression(italic(T[a])),ylab=expression(italic(M)),type="n",xlim=c(0,35),ylim=c(0,1.5),cex.lab=1.5)


  segments(27,0.8,32,0.8,lwd=2,col="coral3")#BMR
  #segments(32,0.8,35,0.95,lwd=1,lty=2)#Thermo_up
  segments(0,2.15,27,0.8,lwd=2,col="coral3")#Therm_down
  segments(2,1,14,0.4,lwd=2, col="cornflowerblue")#Therm_torp
  axis(1,c(2,14,27,31,35),c("",expression(paste("T"["t"])),expression(italic(T[lc])),expression(italic(T[m])),""),cex.axis=1.5)
  axis(2,c(0.2,0.4,0.8,1.5),c("",expression(italic(TMR)),expression(italic(BMR)),""),cex.axis=1.5)

  beta=log(0.8/0.4)/(31-14)
  alpha=0.4/exp(beta*14)
  funtorp <- function(x) alpha*exp(beta*(x))
  points(funtorp(seq(14,27,length.ou=100))~seq(14,27,length.ou=100),type="l",lwd=2, col="cornflowerblue")
  #points(funtorp(seq(27,31,length.ou=100))~seq(27,31,length.ou=100),type="l",lwd=1,lty=2, col="cornflowerblue")


  plot(DATA[,2],results[,1],frame=F)
  abline(0,1)
  abline(10.2184,0.6452)
  summary(lm(scale(DATA[,2])~-1+scale(results[,1])))
  summary(lm(DATA[,2]~results[,1]))
  cor.test(DATA[,2],results[,1],paired = T)
  t.test(DATA[,2],results[,1],paired = T)
  median(DATA[,2]- results[,1],na.rm=T)


  plot(log10(DATA[,4]),log10(results[,2]),frame=F)
  abline(0,1)
  summary(lm(DATA[,4]~-1+results[,2]))
  cor.test(DATA[,4],results[,2],paired = T)
  mean(DATA[,4]- results[,2],na.rm=T)


  cor.test(results[,3],DATA$TaTorpThermo,paired=T)
  table(!is.na(results[,3]),!is.na(DATA$TaTorpThermo))

  success <- rep(NA,28)
  success2 <- rep(NA,28)
  liste2 <- list.files(pattern="*.txt")
  liste2 <- liste2[-c(7:9)]
  valid <-  rep(NA,28)
  diff0 <- 0
  diff12 <- 0
  study <- 0
  N <- rep(NA,28)
  for(t in  1:length(liste2)){
    data <- read.table(liste2[t],header=T)
    data$G <-gsub("1","TOR",data$G)
    data$G <-gsub("2","EUT",data$G)
    data$G <-gsub("3","TNZ",data$G)
    #proportion of valid assignation
    valid[t] <- 1-nrow(data[data$G==0,])/nrow(data)
    #success with only valid data
    success[t] <- table(data[data$G!=0,5]==data[data$G!=0,]$G)["TRUE"]/nrow(data[data$G!=0,])
    #proportion of mismatch due to TLC differences between model and author
    success2[t] <-
      ifelse(
        is.na(table(
          data[data$G!=0&(data$G=="TNZ"|data[,5]=="TNZ"),5]
          !=data[data$G!=0&(data$G=="TNZ"|data[,5]=="TNZ"),]$G)["TRUE"]),
        0,
        table(
          data[data$G!=0&(data$G=="TNZ"|data[,5]=="TNZ"),5]
          !=data[data$G!=0&(data$G=="TNZ"|data[,5]=="TNZ"),]$G)["TRUE"]/nrow(data[data$G!=0,]))

    N[t] <- nrow(data)

    diff0 <- c(diff0,results[t,1]-data[data$G==data$state&data$state!="TNZ"&data$G!="TNZ",2])
    diff12 <- c(diff12,results[t,1]-data[data$G!=data$state&data$state!="TNZ"&data$G!="TNZ",2])
    study <- c(study,rep(t,nrow(data[data$state!="TNZ"&data$G!="TNZ",])))

  }

  #number of studies where no mismatches
  length(success[success==1])
  #number of studies where mismatches are not caused by tlc
  table(success2/(1-success)==0)
  #proportion of mismatches caused by tlc
  success2/(1-success)
  median((success2/(1-success))[(success2/(1-success))!=0],na.rm=T)

  study <- study[-1]
  diff0 <- diff0[-1]
  diff12 <- diff12[-1]
  TAdiff <- c(diff0,diff12)
  Gdiff <- c(rep(0,length(diff0)),rep(1,length(diff12)))
 plot(Gdiff~TAdiff,
      pch="|",
      cex=0.5,
      frame=F,
      ylab="Assign. mismatch",
      xlab=expression(italic(T["lc"])*"-"*italic(T["a"])),
      yaxt="n",
      ylim=c(0,1),
      xlim=c(0,45))
  axis(2,at=c(0,1),labels=c(0,1),lwd=0,las=2)

library(lme4)
summary(mod <-  glmer(Gdiff~TAdiff+(1|as.factor(study)),family="binomial"))
Anova(mod)
expit <- function(x)1/(1+exp(-x))
fun <- function(x)expit(fixef(mod)[1]+fixef(mod)[2]*x)
par(new=T)
plot(fun,
     ylim=c(0,6e-7),
     xlim=c(0,45),
     frame=F,
     ylab="",
     xlab="",
     xaxt="n",
     yaxt="n")
  hist(success)
  mean(success)

  	liste[16]


  	p <- rep(0,1000)

  	for(i in 1:1000){

  	  p[i]  <- binom.test(seq(29000,30000)[i],60000) $ p.value
  	}

  	plot(p[1:1000]~seq(29000,30000)[1:1000],type="l",ylim=c(0,0.1))

  	min(seq(29000,30000)[p>0.05])
  	29760/60000
