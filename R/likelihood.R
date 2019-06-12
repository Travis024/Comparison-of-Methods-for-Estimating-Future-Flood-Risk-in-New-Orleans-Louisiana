## likelihood.R

#===============================================================================
# log(prior) for gev model
#===============================================================================
#
log_prior_gev <- function(parameters,
                          parnames,
                          priors
){
  lpri <- 0
  for (par in parnames) {
    parameter.value <- as.numeric(parameters[match(par,parnames)])
    if(priors[[par]]$type=='normal') {
      lpri <- lpri + dnorm(x=parameter.value, mean=priors[[par]]$mean, sd=priors[[par]]$sd, log=TRUE)
    } else if(priors[[par]]$type=='gamma') {
      lpri <- lpri + dgamma(x=parameter.value, shape=priors[[par]]$shape, rate=priors[[par]]$rate, log=TRUE)
    } else if(priors[[par]]$type=='uniform') {
      lpri <- lpri + dunif(x=parameter.value, min=priors[[par]]$lower, max=priors[[par]]$upper, log=TRUE)
    }
  }
  return(lpri)
}
#===============================================================================


#===============================================================================
# log(likelihood) for gev model
#===============================================================================
log_like_gev <- function(parameters,     # set of GEV parameters (possibly nonstationary)
                         parnames,       # gev_models[[m]]$parnames
                         data_calib,     # processing_output[,"lsl_max"]
                         auxiliary=NULL  # time series of same length as data_calib
){
  if ("mu0" %in% parnames) {
    # location parameter nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    mu <- mu0 + mu1*auxiliary
  } else {
    mu <- parameters[match('mu',parnames)]
  }
  if ("sigma0" %in% parnames) {
    # scale parameter nonstationary
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    sigma <- exp(sigma0 + sigma1*auxiliary)
  } else {
    sigma <- parameters[match('sigma',parnames)]
  }
  if ("xi0" %in% parnames) {
    # shape parameter nonstationary
    xi0 <- parameters[match('xi0',parnames)]
    xi1 <- parameters[match('xi1',parnames)]
    xi <- xi0 + xi1*auxiliary
  } else {
    xi <- parameters[match('xi',parnames)]
  }
  llik <- sum(devd(data_calib, loc=mu, scale=sigma, shape=xi, log=TRUE, type='GEV'))
  return(llik)
}
#===============================================================================


#===============================================================================
# negative log(likelihood) for gev model
#===============================================================================
neg_log_like_gev <- function(parameters,     # set of GEV parameters (possibly nonstationary)
                             parnames,       # gev_models[[m]]$parnames
                             data_calib,     # processing_output[,"lsl_max"]
                             auxiliary=NULL  # time series of same length as data_calib
){
  return(-log_like_gev(parameters,parnames,data_calib,auxiliary))
}
#===============================================================================


#===============================================================================
# log(posterior) for gev model
#===============================================================================
log_post_gev <- function(parameters, parnames, data_calib, priors, auxiliary){
  
  lpost <- 0
  llik <- 0
  lpri <- 0

  # calculate prior
  lpri <- log_prior_gev(parameters=parameters,
                        parnames=parnames,
                        priors=priors)
  
  if(is.finite(lpri)){
    # calculate likelihood (only if parameters pass the prior test)
    llik <- log_like_gev(parameters=parameters,
                         parnames=parnames,
                         data_calib=data_calib,
                         auxiliary=auxiliary)
  }
  
  lpost <- lpri + llik
  return(lpost)
}
#===============================================================================


#===============================================================================
# End
#===============================================================================
