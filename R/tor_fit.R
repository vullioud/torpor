##############################################################################

#'complete_args // internal
#'
#'this function is used internally for tor_fit
#'From Alex Courtiol - informed. He does not want to be listed as an author -
#'
#'@name complete_args
#'@aliases complete_args
#'@param fn a function

.complete_args <- function(fn) {
  ## This function should not be called by the user.
  ## It keeps the default list elements when
  ## a new list with fewer elements is provided.
  env <- parent.frame()
  args <- formals(fn)
  for (arg_name in names(args)) {
    if (is.call(arg <- args[[arg_name]])) {
      if (arg[1] == "list()") {
        arg_input <- mget(names(args), envir = env)[[arg_name]]
        arg_full  <- eval(formals(fn)[[arg_name]])
        if (is.null(names(arg_input))) {
          if (length(arg_input) == length(arg_full)) {
            names(arg_input) <- names(arg_full)
          }
          else {
            stop(paste("The list", arg_name, "should contain names, or be of length equal to the default."))
          }
        }
        arg_full_updated <- utils::modifyList(arg_full, arg_input)
        assign(arg_name, arg_full_updated, envir = env)
      }
    }
  }
  return(NULL)
}


##############################################################################

#'step_1
#'
#' This function estimate the lowest TLC and the bmr
#'
#'@name step_1
#'@param Y A vector of methabolic measure.
#'@param Ta A vector of ambient temperature.
#'@export
#'@examples
#'step_1(test_data$VO2, test_data$Ta)

step_1 <- function(Y,Ta){


  SEQX <- seq(from = max(Ta,na.rm = TRUE), to = min(Ta,na.rm = TRUE), length.out = 100)
  Nb <- rep(NA, 100)
  whitetest <- rep(NA, 100)
  p <- rep(NA, 100)

  #Define heterosc and slope
  for (j in 1:100) {

      Nb[j] <- length(Y[Ta > SEQX[j]])
    if (Nb[j] < 11) {
      whitetest[j] <- NA
    p[j] <- NA
    } else {
      whitetest[j] <- lmtest::bptest(stats::lm(Y[Ta > SEQX[j]] ~ Ta[Ta > SEQX[j]]))$p.value

      mod <- stats::lm(Y[Ta> SEQX[j]]~Ta[Ta> SEQX[j]])
      p[j] <- ifelse(stats::coefficients(mod)[2] > 0,NA, summary(mod)$coefficients[2,4])}
  }

  #Define "low TLC" and first BMR
  #if (length(SEQX[p < 0.01 & !is.na(p)]) == 0) stop("Error: not possible to estimate MTNZ or TLC")


  test_1 <- ifelse(length(stats::na.omit(SEQX[whitetest > 0.05])) >= 1, min(SEQX[whitetest > 0.05], na.rm = TRUE), NA)
  test_2 <- ifelse(length(stats::na.omit(SEQX[p < 0.01])) >= 1, SEQX[match(max(SEQX[p < 0.01],na.rm = TRUE) ,SEQX) - 1], NA)

  heterosc <- ifelse(length(stats::na.omit(c(test_1, test_2))) < 1, NA, max(test_1, test_2, na.rm = TRUE))

  if (is.na(heterosc)) stop("TLC and BMR not estimable, please provide them by hand")

  if (length(Y[Ta > heterosc]) < 10) warning("MTNZ computed on less than 10 points")
  bmr <- mean(Y[Ta > heterosc], na.rm = TRUE)

  out <- c(bmr,heterosc)


names(out) <- c("bmr", "low_tlc")
return(out)
}

#'step_1_and_2
#'
#' This function estimate the lowest TLC and the bmr
#'
#'@name step_1_and_2
#'@param Y A vector of methabolic measure.
#'@param Ta A vector of ambient temperature.
#'@inheritParams tor_fit
#'@export
step_1_and_2 <- function(Y, Ta, fitting_options = list(ni = 500000,
                                                       nt = 10,
                                                       nb = 300000,
                                                       nc = 3)) {

  out_1 <- step_1(Y, Ta) ## run first step

  bmr <- out_1["bmr"] ## extract bmr
  heterosc <- out_1["low_tlc"] ## extract low tlc

  Ym <- mean(Y, na.rm = TRUE) ## keep the mean
  Ta2 <- Ta[Ta < heterosc]
  Y2 <- Y[Ta < heterosc]

  inits <- list(
      tauy = runif(1,0.1,1),
      coef = runif(1,1,1.1),
      p = runif(2),
      Tbe = runif(1,heterosc,50),
      tlc = runif(1,heterosc,max(Ta)),
      Tt = runif(1),
      TMR = runif(1,0,bmr*0.8/Ym),
      intc1 = runif(1,0,bmr*0.8/Ym))

  inits_hetero_list <- rep(list(inits), fitting_options[["nc"]])

  params_hetero <- c("tlc","betat","inte")

  win_data <- list(Y = Y2/Ym,
                   NbObservations = length(Y2),
                   Ta = Ta2,
                   heterosc = heterosc,
                   BMR = bmr/Ym,
                   Max = max(Ta),
                   Ym = Ym)

  # MCMC settings


  path_to_model_1 <- system.file("extdata", "hetero_1.txt",  package = "torpor")

  out <- jagsUI::jags(data = win_data,
                        inits = inits_hetero_list,
                        parameters.to.save = params_hetero,
                        model.file = path_to_model_1,
                        n.chains = fitting_options[["nc"]],
                        n.thin = fitting_options[["nt"]],
                        n.iter = fitting_options[["ni"]],
                        n.burnin = fitting_options[["nb"]],
                        parallel = TRUE,
                        verbose = FALSE,
                        store.data = TRUE)

return(out)
}

#'step_1_and_2
#'
#' This function estimate the lowest TLC and the bmr
#'
#'@name step_3
#'@export
#'@param bmr bmr value if not estimated
#'@param tlc tlc value if not estimated
#'@inheritParams step_1_and_2
#'@examples
#'\dontrun{
#'step_3(Ta = test_data$Ta, Y = test_data$VO2)
#'}
step_3 <- function(Y, Ta, bmr = NULL, tlc = NULL, fitting_options = list(ni = 500000,
                                                                         nt = 10,
                                                                         nb = 300000,
                                                                         nc = 3)){

  .complete_args(step_3)

  ## check input
  if (length(Y) != length(Ta)) stop("Ta and MR have not the same length")

  data <- cbind(Y, Ta)[!is.na(Y) & !is.na(Ta), ] ## remove Nas

  Y <- data[, "Y"]
  Ta <- data[, "Ta"]

  ## run step 1 and 2 to get bmr and tlc

  if (is.null(bmr) | is.null(tlc)) {
    message("BMR and TLC are being
            estimated from the data")
    out_2 <- step_1_and_2(Ta = Ta, Y = Y)
    bmr <- mean(Y[Ta > out_2$mean$tlc])*out_2$data$Ym
    tlc <- out_2$mean$tlc*out_2$data$Ym
  }

  Ym <- mean(Y, na.rm = T)
  heterosc <- out_2$data$heterosc
  ## initial values
  inits_2 <- list(
    tauy = runif(1,0.1,1),
    coef = runif(1,1,1.1),
    p = runif(2),
    Tbe = runif(1,heterosc,50),
    tlc = runif(1,heterosc,max(Ta)),
    Tt = runif(1),
    TMR = runif(1,0,bmr*0.8/Ym),
    intc1 = runif(1,0,bmr*0.8/Ym))

  inits_hetero_list_2 <- rep(list(inits_2), fitting_options[["nc"]])

  params_hetero_2 <- c("tau",
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

  win.data_2 <- list(Y = Y/Ym,
                     NbObservations = length(Y),
                     Ta = Ta,
                     BMR = bmr/Ym,
                     tlc = tlc,
                     Ym = Ym)

  path_to_model_2 <- system.file("extdata", "hetero_2.txt",  package = "torpor")

  out2 <- jagsUI::jags(data = win.data_2,
                       inits = inits_hetero_list_2,
                       parameters.to.save = params_hetero_2,
                       model.file = path_to_model_2,
                       n.chains = fitting_options[["nc"]],
                       n.thin = fitting_options[["nt"]],
                       n.iter = fitting_options[["ni"]],
                       n.burnin = fitting_options[["nb"]],
                       parallel = TRUE,
                       verbose = FALSE,
                       store.data = TRUE)

  return(out2)

}


#'step_3_bis
#'
#' This function estimate the lowest TLC and the bmr
#'
#'@name step_3_bis
#'@param mod A fitted model
#'@export
step_3_bis <- function(mod){
  Ta <- mod$data$Ta
  Y <- mod$data$Y

  G <- mod$mean$G
  limit <- .get_limit_value(nrow(mod$samples[[1]]))
  G[sqrt((G - round(G))^2) > limit] <- 0
  G <- round(G)

  if (length(G[G == 1]) != 0 ) {

    for (i in 1:length(Y <= tlc)) {

      #Only proceed to this step if there are at least one torpid (G==1) value.
      G[i] <- ifelse(Y[i] <= stats::median(tor_predict_fun(x = Ta[i],
                                                           Tt = mod$sims.list$Tt,
                                                           intr = mod$sims.list$intr,
                                                           intc = mod$sims.list$intc,
                                                           betat = mod$sims.list$betat,
                                                           betac = mod$sims.list$betac,
                                                           Ym = mod$data$Ym)) & G[i] != 3, 1, G[i])
      G[i] <- ifelse(Y[i] >=  stats::median(eut_predict_fun(x = Ta[i],
                                                            inte = mod$sims.list$inte,
                                                            betat = mod$sims.list$betat,
                                                            Ym = mod$data$Ym)) & G[i] != 3, 2 ,G[i])
    }

  }
  return(G)
}

#' step_4
#'
#'[step_4()] fits a binomial mixture model using Bayesian
#'inference.  The function considers the assumed relation between metabolic rate (MR) and ambient temperature (Ta)
#'(Speakman & Thomas 2003). In the hypothermic state (torpor) and above some
#'threshold Ta (Tmin), MR follows an exponential curve reflecting the Arrhenius
#'rate enhancing effect of temperature on chemical reactions, whereas below
#'Tmin, it increases linearly with decreasing Ta to maintain a minimal body temperature (Tb) in
#'torpor. In the euthermic state, MR solely increases linearly with decreasing Ta.
#'The function uses Rjags in the background and enables users to specify
#'some - but not all - sampling parameters. The structure of the model can also
#'be changed. Users who want more flexibility are encouraged to use Rjags
#'directly.
#'More information about the model can be found in the \code{vignette("model_description", package = "torpor")}.
#'
#'This model should be applied with sufficient sample size and if
#'evidence suggest that individuals under study will conform to the previously
#'described pattern while in torpor.
#'
#'@name step_4
#'@inheritParams step_3
#'@export
#'@examples
#'\dontrun{
#'test_mod <- step_4(Ta = test_data3$Ta, Y = test_data3$VO2)
#'}
step_4 <- function(Ta, Y,  bmr = NULL, tlc = NULL, fitting_options = list(ni = 5000,
                                                                          nt = 10,
                                                                          nb = 3000,
                                                                          nc = 3), .debug = FALSE) {


  if (length(Ta) != length(Y)) stop("Y, Ta don´t have the same length")
  ## run step 1-3
  out_2 <- step_3(Ta = Ta, Y = Y,  bmr = bmr, tlc = tlc, fitting_options = fitting_options)


  ## manual
  G <- step_3_bis(out_2)
  if (.debug) return(list(out_2, G))

  Y <- out_2$data$Y*out_2$data$Ym
  Ta <- out_2$data$Ta
  tlc <- out_2$data$tlc*out_2$data$Ym
  bmr <- out_2$data$BMR*out_2$data$Ym
  Ym <- out_2$data$Ym
  win.data_3 <- list(Y = Y[G != 0 & G != 3] / Ym,
                     NbObservations = length(Y[G != 0 & G != 3]),
                     Ta = Ta[G != 0 & G != 3],
                     BMR = bmr / Ym,
                     tlc = tlc,
                     G = G[G != 0 & G != 3],
                     Ym = Ym)


  path_to_model_3 <- system.file("extdata", "hetero_3.txt",  package = "torpor")

  inits_3 <- list(
    tauy = runif(1,0.1,1),
    coef = runif(1,1,1.1),
    p = runif(3),
    Tbe = runif(1,tlc,50),
    intc1 = runif(1,0,bmr*0.8/Ym),
    Tt = runif(1),
    TMR = runif(1,0,bmr*0.8/Ym))

  inits_hetero_list_3 <- rep(list(inits_3), fitting_options[["nc"]])

  out3 <- jagsUI::jags(data = win.data_3,
                       inits =  inits_hetero_list_3,
                       parameters.to.save = out_2$parameters,
                       model.file = path_to_model_3,
                       n.chains = fitting_options[["nc"]],
                       n.thin = fitting_options[["nt"]],
                       n.iter = fitting_options[["ni"]],
                       n.burnin = fitting_options[["nb"]],
                       parallel = TRUE,
                       verbose = FALSE,
                       store.data = TRUE)


  return(out3)
}


#'.get_limit_value
#'
#' find the limit value
#'
#'@name .get_limit_value
#'@param nb_samples size of the samples
.get_limit_value <- function(nb_samples){

  Seq <- round(nb_samples*0.9/2):round(nb_samples/2)

  Pvalues <- rep(NA,length(Seq))



  for(i in 1:length(Seq)){

    Pvalues[i] <- binom.test(Seq[i],nb_samples)$p.value}


  max(Seq[Pvalues<0.05])/nb_samples
}
