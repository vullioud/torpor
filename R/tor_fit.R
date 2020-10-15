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

#'Find the low Tlc.
#'
#' This function estimate the lowest Tlc and the mtnz
#'
#'@name find_low_tlc_mtnz
#'@param M A vector of methabolic measure.
#'@param Ta A vector of ambient temperature.
#'@return A named vector with mtmz and low_tlc.
#'@export
#'@examples
#'t_1 <- find_low_tlc_mtnz(M = test_data$VO2, Ta = test_data$Ta)
find_low_tlc_mtnz <- function(M,Ta){
  ### clean arguments function
  Y <- M
  x <- stats::na.omit(cbind(Ta, Y))

  Y <- x[,2]
  Ta <- x[,1]

  SEQX <- seq(from = max(Ta, na.rm = TRUE), to = min(Ta, na.rm = TRUE), length.out = 100)

  Nb <- rep(NA, 100)
  whitetest <- rep(NA, 100)
  p <- rep(NA, 100)
  p2 <- rep(NA, 100)
  #Define low_tlc and slope
  for (j in 1:100) {

      Nb[j] <- length(Y[Ta > SEQX[j]])

    if (Nb[j] < 11) {

    whitetest[j] <- NA
    p[j] <- NA
    p2[j] <- NA

    } else {
      whitetest[j] <- lmtest::bptest(stats::lm(Y[Ta > SEQX[j]] ~ Ta[Ta > SEQX[j]]))$p.value

      mod <- stats::lm(Y[Ta > SEQX[j]] ~ Ta[Ta> SEQX[j]])
      p[j] <- ifelse(stats::coefficients(mod)[2] > 0, NA, summary(mod)$coefficients[2,4])
      p2[j] <- ifelse(stats::coefficients(mod)[2] < 0, NA, summary(mod)$coefficients[2,4])


      }
  }

  if (length(SEQX[p2 < 0.01 & !is.na(p2)]) > 0 &
      length(SEQX[p < 0.01 & !is.na(p)]) == 0 &
      length(SEQX[whitetest < 0.01]) == 0) {

    stop("Tlc and Mtnz can not be estimated: provide values")

  } else {

  test_1 <- ifelse(length(stats::na.omit(SEQX[whitetest > 0.05])) >= 1, min(SEQX[whitetest > 0.05], na.rm = TRUE), NA)

  test_2 <- ifelse(length(stats::na.omit(SEQX[p < 0.01])) >= 1, SEQX[match(max(SEQX[p < 0.01],na.rm = TRUE) ,SEQX) - 1], NA)

  low_tlc <- ifelse(length(stats::na.omit(c(test_1, test_2))) == 0, NA, max(test_1, test_2, na.rm = TRUE)) ### check if ifelse necessary

  }
  if (is.na(low_tlc)) stop("Tlc and Mtnz can not be estimated: provide values")
  if (length(Y[Ta > low_tlc]) < 10) warning("Mtnz computed on less than 10 points")

  MTNZ <- mean(Y[Ta > low_tlc], na.rm = TRUE)
  out <- c(MTNZ,low_tlc)
  names(out) <- c("MTNZ", "low_tlc")
  out
}

#' Estimate Tlc
#'
#' This function estimate the lowest Tlc and the MTNZ
#'
#'@name estimate_tlc_mtnz
#'@param M A vector of methabolic measure.
#'@param Ta A vector of ambient temperature.
#'@param MTNZ A value of MTNZ.
#'@param fitting_options a list of fitting option to pass to jags.
#'@export
#'@examples
#'t <- estimate_tlc_mtnz(M = test_data2$VO2ms, Ta = test_data2$Ta, MTNZ = NULL)

estimate_tlc_mtnz <- function(M, Ta, MTNZ, fitting_options = list(ni = 50000,
                                                       nt = 10,
                                                       nb = 20000,
                                                       nc = 2,
                                                       parallel = TRUE)) { ## add parallel to fitting options
  .complete_args(estimate_tlc_mtnz)
  Y <- M
  MTNZ_input <- MTNZ


  out_1 <- find_low_tlc_mtnz(Y, Ta) ## run first step

  if(!is.null(MTNZ_input)) {
    MTNZ <- MTNZ_input
  } else {
  MTNZ <- out_1["MTNZ"]
  } ## extract MTNZ
  low_tlc <- out_1["low_tlc"] ## extract low tlc


  Ym <- mean(Y, na.rm = TRUE) ## keep the mean
  Ta2 <- Ta[Ta < low_tlc]
  Y2 <- Y[Ta < low_tlc]

  inits <- list(
      tauy2= stats::runif(1,0.05,0.1),
      tauy1= stats::runif(1,0.03,0.05),
      p = stats::runif(2),
      Tbe = stats::runif(1,low_tlc,50),
      Tlc = stats::runif(1, low_tlc, max(Ta)),
      Tbt = stats::runif(1, 0, max(Ta2)),
      TMR = stats::runif(1,1e-5, MTNZ*0.8/Ym),
      MR = stats::runif(1, MTNZ*0.8/Ym, MTNZ/Ym))

  inits_hetero_list <- rep(list(inits), fitting_options[["nc"]])

  params_hetero <- c("Tlc","betat","inte")

  win_data <- list(Y = Y2/Ym,
                   NbObservations = length(Y2),
                   Ta = Ta2,
                   MTNZ = MTNZ/Ym,
                   Max = max(Ta))

  # MCMC settings
  path_to_model_1 <- system.file("extdata", "hetero1.txt",  package = "torpor")

  mod <- jagsUI::jags(data = win_data,
                        inits = inits_hetero_list,
                        parameters.to.save = params_hetero,
                        model.file = path_to_model_1,
                        n.chains = fitting_options[["nc"]],
                        n.thin = fitting_options[["nt"]],
                        n.iter = fitting_options[["ni"]],
                        n.burnin = fitting_options[["nb"]],
                        parallel = fitting_options[["parallel"]],
                        verbose = FALSE,
                        store.data = TRUE)

  ## extract the important output
  tlc_estimated <- stats::median(mod$sims.list$Tlc) ## estimated tlc
  tlc_distribution <- mod$sims.list$Tlc ## distribution of tlc

  if(is.null(MTNZ_input)) {
  mtnz_estimated <- mean(Y[Ta >= tlc_estimated], na.rm = TRUE) # estimated MTNZ (mean of the points)
  mtnz_points <- stats::na.omit(Y[Ta >= tlc_estimated])
  } else {
    mtnz_estimated <- MTNZ_input
    mtnz_points <- MTNZ_input
  }

if(length(stats::na.omit(Y[Ta >= tlc_estimated])) < 10 & is.null(MTNZ_input)) warning("Mtnz computed on less than 10 points")

return(list(model_1 = mod,
            tlc_estimated = tlc_estimated,
            tlc_distribution = tlc_distribution,
            mtnz_estimated = mtnz_estimated,
            mtnz_points = mtnz_points,
            Ta2 = Ta2))
}

#' Estimate assignation
#'
#' This function estimate the assignation into torpor or euthermy. It is used
#' internally prior to estimate the parameters
#'
#'@name estimate_assignation
#'@export
#'@inheritParams estimate_tlc_mtnz
#'@param MTNZ mtnz value if not estimated
#'@param Tlc tlc value if not estimated
#'@examples
#'\dontrun{
#'t2_no_input <- estimate_assignation(Ta = test_data2$Ta, M = test_data2$VO2)
#'to_input <- estimate_assignation(Ta = test_data2$Ta, M = test_data2$VO2, MTNZ = 1.49, Tlc = 28.8)
#'}
estimate_assignation <- function(M, Ta,
                                 MTNZ = NULL, Tlc = NULL,
                                 fitting_options = list(ni = 50000,
                                                        nt = 10,
                                                        nb = 20000,
                                                        nc = 2,
                                                        parallel = TRUE)){

  runif <- NULL
  .complete_args(estimate_assignation)

  ## check input
  Y <- M
  if (length(Y) != length(Ta)) stop("Ta and M have not the same length")

  data <- cbind(Y, Ta)[!is.na(Y) & !is.na(Ta), ] ## remove Nas

  Y <- data[, "Y"]
  Ta <- data[, "Ta"]
  Ym <- mean(Y, na.rm = T)

  ## run step 1 and 2 to get mtnz and Tlc

  if (is.null(MTNZ) & (is.null(Tlc))) { ### need to change the filter it runs all the time

    message("Mtnz and Tlc are being estimated from the data")

  out_mtnz_tlc <- estimate_tlc_mtnz(Ta = Ta, M = Y, MTNZ = MTNZ, fitting_options = fitting_options)

  } else if (!(is.null(MTNZ)) & (is.null(Tlc))) { ### need to change the filter it runs all the time

    message("Tlc is being estimated from the data")

    out_mtnz_tlc <- estimate_tlc_mtnz(Ta = Ta, M = Y, MTNZ = MTNZ, fitting_options = fitting_options)

  } else { ## keep the structure of output of estimate_tlc_mtnz when there are provided as input

    out_mtnz_tlc <- list(model_1 = NULL,
         tlc_estimated = Tlc,
         tlc_distribution = Tlc,
         mtnz_estimated = MTNZ,
         mtnz_points = MTNZ,
         Ta2 = Ta[Ta < Tlc])
  }
  Tlc <- out_mtnz_tlc$tlc_estimated
  MTNZ <- out_mtnz_tlc$mtnz_estimated

  ## initial values
  inits_2 <- list(
    tauy2=stats::runif(1,0.05,0.1),
    tauy1=stats::runif(1,0.03,0.05),
    p = stats::runif(3),
    Tbe = stats::runif(1,Tlc, 50),
    Tbt = stats::runif(1, Tlc - 1, Tlc),
    TMR = stats::runif(1,1e-5,MTNZ*0.8/Ym),
    MR = stats::runif(1, MTNZ*0.8/Ym, MTNZ/Ym))

  inits_hetero_list_2 <- rep(list(inits_2), fitting_options[["nc"]])

  params_hetero_2 <- c("tau","G","p","inte","intc","intr","betat","betac","Tt","TMR","MR","tauy1", "tauy2")

  win.data_2 <- list(Y = Y/Ym,
                     NbObservations = length(Y),
                     Ta = Ta,
                     MTNZ = MTNZ/Ym,
                     Tlc = Tlc)

  path_to_model_2 <- system.file("extdata", "hetero2.txt",  package = "torpor")

  mod_assignation <- jagsUI::jags(data = win.data_2,
                       inits = inits_hetero_list_2,
                       parameters.to.save = params_hetero_2,
                       model.file = path_to_model_2,
                       n.chains = fitting_options[["nc"]],
                       n.thin = fitting_options[["nt"]],
                       n.iter = fitting_options[["ni"]],
                       n.burnin = fitting_options[["nb"]],
                       parallel = fitting_options[["parallel"]],
                       verbose = FALSE,
                       store.data = TRUE)

  out_assignation <- list(
       mod_assignation = mod_assignation,
       out_mtnz_tlc = out_mtnz_tlc,
       data = list(Y = Y,
                   Ta = Ta,
                   Ym = Ym))

  out_assignation$assignation <- .step_3_bis(out_assignation = out_assignation,
                             fitting_options = fitting_options)

  out_assignation
}


#'step_3_bis
#'
#' This function estimate the lowest TLC and the MTNZ
#'
#'@name .step_3_bis
#'@param out_assignation output of estimatate_assignation.
#'@param fitting_options fitting options
.step_3_bis <- function(out_assignation, fitting_options){
 ## comes from hetero2
  Ta <- out_assignation$data$Ta
  Y <- out_assignation$data$Y

  nc <- fitting_options[["nc"]]
  nt <- fitting_options[["nt"]]
  ni <- fitting_options[["ni"]]
  nb <- fitting_options[["nb"]]
  ## extract model assignation
  mod_assignation <- out_assignation$mod_assignation

  G <- mod_assignation$mean$G ## extract the G


  nbsamples <- ((ni - nb)*nc)/nt ## look how to take them from mod2

  ## inside small fun
  expit <- function(x){1/(1 + exp(-x))} ##step not to be done if Tlc and MTNZ and given
  funabove <- function(x, mod){expit(stats::coefficients(mod)[1] + stats::coefficients(mod)[2]*x)} ##step not to be done if Tlc and MTNZ and given
  funbelow <- function(x, mod){-(funabove(x, mod)- 1)} ##step not to be done if Tlc and MTNZ and given

  ####
  probStatus <- 1 - sqrt((G - round(G))^2)

  length_tlc_dist <- length(out_assignation$out_mtnz_tlc$tlc_distribution)

  if(length_tlc_dist > 1) {
  mod_glm <- stats::glm(stats::pnorm(q = out_assignation$out_mtnz_tlc$tlc_distribution, ## changed for Tlc distribution
                       mean = out_assignation$out_mtnz_tlc$tlc_estimated,
                       sd = ifelse(length_tlc_dist > 1, stats::sd(out_assignation$out_mtnz_tlc$tlc_distribution), 0)) ~
                   out_assignation$out_mtnz_tlc$tlc_distribution,
                 family = "binomial") ##step not to be done if Tlc and MTNZ and given


  probStatus <- probStatus*ifelse(G == 3, funabove(Ta, mod_glm), funbelow(Ta, mod_glm))  ##step not to be done if Tlc and MTNZ and given
  }
  pstatus <- rep(NA,length(Ta))

  for(i in 1:length(Ta)) {

    pstatus[i] <- as.numeric(stats::binom.test(round(nbsamples*probStatus[i]),
                                        nbsamples,
                                        alternative = "greater", p = 0.5)$p.value)

  }

  G[pstatus>0.01] <- 0
  G <- round(G)

  if(length(G[G ==1])!=0){

    for(i in 1:length(G)){
    #Only proceed to this step if there are at least one torpid (G==1) value.
    G[i] <- ifelse(Y[i] <= stats::median(tor_predict_fun(x = Ta[i],
                                                         Tt = mod_assignation$sims.list$Tt,
                                                         intr = mod_assignation$sims.list$intr,
                                                         intc = mod_assignation$sims.list$intc,
                                                         betat = mod_assignation$sims.list$betat,
                                                         betac = mod_assignation$sims.list$betac,
                                                         Ym = out_assignation$data$Ym)) & G[i] != 3, 1, G[i])

    G[i] <- ifelse(Y[i] >=  stats::median(eut_predict_fun(x = Ta[i],
                                                          inte = mod_assignation$sims.list$inte,
                                                          betat = mod_assignation$sims.list$betat,
                                                          Ym = out_assignation$data$Ym)) & G[i] != 3, 2 ,G[i])
    }
  }

  list(G = G,
       pstatus = pstatus)

}


#' Fit a binomial mixture model
#'
#'[tor_fit()] fits a binomial mixture model using Bayesian
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
#'More information about the model can be found in Fassel et al. (in prep).
#'
#'This model should be applied with sufficient sample size and if
#'evidence suggest that individuals under study will conform to the previously
#'described pattern while in torpor.
#'
#'@name tor_fit
#'@inheritParams estimate_assignation
#'@return A list of class tor_obj
#'@export
#'@examples
#'\dontrun{
#'test_mod <- tor_fit(Ta = test_data3$Ta, M = test_data3$Y,
#'                    fitting_options = list(parallel = TRUE))
#'test_mod2 <- tor_fit(Ta = test_data2$Ta, M = test_data2$VO2ms, MTNZ = 1.8, Tlc = 29.2,
#'                    fitting_options = list(parallel = TRUE))
#'test_mod3 <- tor_fit(Ta = test_data2$Ta, M = test_data2$VO2ms, MTNZ = 1.8,
#'                    fitting_options = list(parallel = TRUE))
#'}
tor_fit <- function(Ta, M,
                   MTNZ = NULL,
                   Tlc = NULL,
                   fitting_options = list(ni = 50000,
                                          nt = 10,
                                          nb = 30000,
                                          nc = 3,
                                          parallel = TRUE)) {


  if (length(Ta) != length(M)) stop("M, Ta do not have the same length")
  ## run step 1-3
  .complete_args(tor_fit)


  out_assignation <- estimate_assignation(M = M, Ta = Ta, MTNZ = MTNZ, Tlc = Tlc, fitting_options)

  G <- out_assignation$assignation$G
  Y <- out_assignation$data$Y
  Ta <- out_assignation$data$Ta
  Tlc <- out_assignation$out_mtnz_tlc$tlc_estimated ## add option to include it manualy.
  MTNZ <- out_assignation$out_mtnz_tlc$mtnz_estimated ## add option to include it manualy.
  Ym <- out_assignation$data$Ym

  win.data_3 <- list(Y = Y[G != 0 & G != 3] / Ym,
                     NbObservations = length(Y[G != 0 & G != 3]),
                     Ta = Ta[G != 0 & G != 3],
                     MTNZ = MTNZ / Ym,
                     Tlc = Tlc,
                     G = G[G != 0 & G != 3])

  path_to_model_3 <- system.file("extdata", "hetero3.txt",  package = "torpor")

  inits_3 <- list(Tbe = stats::runif(1,Tlc, 50),
                  Tbt = stats::runif(1, Tlc - 1, Tlc),
                  TMR = stats::runif(1,1e-5,MTNZ*0.8/Ym),
                  MR = stats::runif(1, MTNZ*0.8/Ym, MTNZ/Ym),
                  tauy2 = stats::runif(1,0.05,0.1),
                  tauy1 = stats::runif(1,0.03,0.05))

  inits_hetero_list_3 <- rep(list(inits_3), fitting_options[["nc"]])
  params <- params_hetero_2 <- c("tau","inte","intc","intr","betat",
                                 "betac","Tt","TMR","MR","tauy1", "tauy2", "Tbt", "Tbe")

  out_4 <- jagsUI::jags(data = win.data_3,
                       inits =  inits_hetero_list_3,
                       parameters.to.save = params,
                       model.file = path_to_model_3,
                       n.chains = fitting_options[["nc"]],
                       n.thin = fitting_options[["nt"]],
                       n.iter = fitting_options[["ni"]],
                       n.burnin = fitting_options[["nb"]],
                       parallel = fitting_options[["parallel"]],
                       verbose = FALSE,
                       store.data = TRUE)

  ### check convergence.... if Rhat > 1.1 warrning sur parameters.
  for(i in seq_along(out_4$Rhat)){
   if (any(out_4$Rhat[i][[1]] > 1.1, na.rm = TRUE)) warning(paste("Rhat > 1.1 on parameter:", names(out_4$Rhat[i])))
  }

  out_assignation$mod_parameter <- out_4

  class(out_assignation) <- c("tor_obj", "list")
  print(out_assignation)
  return(invisible(out_assignation))

}

### coommon method to ease the usage
#######
#' @export
#' @method print tor_obj
print.tor_obj <- function(x, ...) {
  print(tor_summarise(x, ...))
}

#######
#' @export
#' @method summary tor_obj
summary.tor_obj <- function(object, ...) {
  tor_summarise(object)
}

#######
#' @export
#' @method plot tor_obj
plot.tor_obj <- function(x, ...) {
  tor_plot(tor_obj = x,...)
}
class(plot.tor_obj)
