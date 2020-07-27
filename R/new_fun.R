Define_uncertain <- function(mod, mod2,  ni, nb, nc, nt){

  ## comes from hetero2
  G <- mod2$mean$G
  Ta <- mod2$data$Ta
  Y <- mod2$data$Y*mod2$data$Ym

  nbsamples <- ((ni - nb)*nc)/nt ## look how to take them from mod2

  ## tlc from mod1
  mod_glm <- glm(pnorm(mod$sims.list$tlc,
                    mod$mean$tlc,
                    sd(mod$sims.list$tlc))~ mod$sims.list$tlc, family="binomial")


  expit <- function(x){1/(1+exp(-x))}

  funabove <- function(x, mod){expit(coefficients(mod)[1]+ coefficients(mod)[2]*x)}

  funbelow <- function(x, mod){-(funabove(x, mod)-1)}

  probStatus <- 1-sqrt((G-round(G))^2)
  probStatus <- probStatus*ifelse(G==3,funabove(Ta, mod_glm),funbelow(Ta, mod_glm))

  pstatus <- rep(NA,length(Ta))

  for(i in 1:length(Ta)){
    pstatus[i] <- as.numeric(binom.test(round(nbsamples*probStatus[i]),
                                        nbsamples,
                                        alternative="greater",p=0.5)$p.value)

  }

  G[pstatus>0.01] <- 0
  G <- round(G)

  if(length(G[G==1])!=0){

  #Only proceed to this step if there are at least one torpid (G==1) value.
  G[i] <- ifelse(Y[i] <= stats::median(tor_predict_fun(x = Ta[i],
                                                       Tt = mod2$sims.list$Tt,
                                                       intr = mod2$sims.list$intr,
                                                       intc = mod2$sims.list$intc,
                                                       betat = mod2$sims.list$betat,
                                                       betac = mod2$sims.list$betac,
                                                       Ym = mod2$data$Ym)) & G[i] != 3, 1, G[i])

  G[i] <- ifelse(Y[i] >=  stats::median(eut_predict_fun(x = Ta[i],
                                                        inte = mod2$sims.list$inte,
                                                        betat = mod2$sims.list$betat,
                                                        Ym = mod2$data$Ym)) & G[i] != 3, 2 ,G[i])
  }
G
}

#Any Y value higher than Ymeann should be considered as Euthermic and Y value lower than Ymeant should be considered as torpid

# for(i in 1:length(Y)){
#
#
#
#   if(length(G[G==1])!=0){
#
#     G[i] <-
#
#       ifelse(Y[i]>= mean(funnorm(Ta[i]))&G[i]!=3,2,G[i])
#
#     G[i] <-
#
#       ifelse(Y[i]<= mean(funtorp(Ta[i]))&G[i]!=3,1,G[i])}
#
# }
