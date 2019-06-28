#' get_overlap
#'
#' Compare the distribution overlap between a prior and a posterior // used internally only
#'@name get_overlap
#'@aliases get_overlap
#'@param mod a fitted model with fit_torpor
#'@param params a parameter of the model
#'@param priors a posterior distribution (as. numerical vector)
#'@importFrom magrittr %>%
#'@return a numeric value
get_overlap <- function(mod, params, priors) {

    chains <- purrr::map_dfr(mod$samples, ~ as.data.frame(.x) %>%
                                dplyr::select(paste(params)))

    as.numeric(round(overlapping::overlap(x = list(chains[,1], priors))$OV, digits = 3))
}

#' tor_overlap
#'
#'The function tor_overlap() generates prior posterior overlap values for Tm,
#'Tt and Betae.
#'Note: Values larger than 0.3 should lead to the conclusion that conforming
#'torpor, regulated torpor or thermoregulation respectively could not be modeled
#' with the data provided (Fasel et al. (in prep)).
#'@name tor_overlap
#'@aliases tor_overlap
#'@param mod a fitted model with fit_torpor
#'@return a list of 3. With overlapping values for tlc, Tt and Betat
#'@export
#'@examples
#'data(test_data2)
#'test <- tor_fit(MR = test_data2[,2],
#'Ta = test_data2[, 1],
#'BMR = 1.49,
#'TLC = 28.8,
#'model = NULL,
#'fitting_options = list(nc = 1, ni = 5000, nb = 3000))
#'tor_overlap(mod = test)

tor_overlap <- function(mod){

prior_tlc <- truncnorm::rtruncnorm(20000,a= mod$mean$TLC,mean=0,sd=100)
prior_Tt <- truncnorm::rtruncnorm(20000,b=mod$mean$TLC,mean=0,sd=100)
prior_betat <- truncnorm::rtruncnorm(20000,b=0,mean=0,sd=100)

tlc_overlap <- get_overlap(mod, "tlc", prior_tlc)
Tt_overlap <- get_overlap(mod, "Tt", prior_Tt)
betat_overlap <- get_overlap(mod, "betat", prior_betat)

if(Tt_overlap >= 0.3){
  warning("Parameters Tt and intr not identifiable")
}
if(tlc_overlap >= 0.3){
  warning("Parameters tlc, betac, and intc not identifiable")
}
if(betat_overlap >= 0.3){
  warning("Parameters tlc, betac, and intc not identifiable")
}

out <- list(tlc = tlc_overlap,
                    Tt = Tt_overlap,
                    Betat =betat_overlap)

return(out)
}

