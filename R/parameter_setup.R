## parameter_setup.R

nmodel <- 8
gev_models <- vector("list", nmodel)

# set up parameter names for each model, and
# set up parameter bounds for each model (need for optimization for priors)

# all stationary
gev_models[[1]]$parnames <- c("mu", "sigma", "xi")
gev_models[[1]]$bound_lower <- c(0, 0, -3)
gev_models[[1]]$bound_upper <- c(5000, 1000, 3)

# mu nonstationary
gev_models[[2]]$parnames <- c("mu0", "mu1", "sigma", "xi")
gev_models[[2]]$bound_lower <- c(0, -500, 0, -3)
gev_models[[2]]$bound_upper <- c(5000, 500, 1000, 3)

# sigma nonstationary
gev_models[[3]]$parnames <- c("mu", "sigma0", "sigma1", "xi")
gev_models[[3]]$bound_lower <- c(0, 0, -200, -3)
gev_models[[3]]$bound_upper <- c(5000, 1000, 200, 3)

# xi nonstationary
gev_models[[4]]$parnames <- c("mu", "sigma", "xi0", "xi1")
gev_models[[4]]$bound_lower <- c(0, 0, -3, -3)
gev_models[[4]]$bound_upper <- c(5000, 1000, 3, 3)

# mu and sigma nonstationary
gev_models[[5]]$parnames <- c("mu0", "mu1", "sigma0", "sigma1", "xi")
gev_models[[5]]$bound_lower <- c(0, -500, 0, -200, -3)
gev_models[[5]]$bound_upper <- c(5000, 500, 1000, 200, 3)

# mu and xi nonstationary
gev_models[[6]]$parnames <- c("mu0", "mu1", "sigma", "xi0", "xi1")
gev_models[[6]]$bound_lower <- c(0, -500, 0, -3, -3)
gev_models[[6]]$bound_upper <- c(5000, 500, 1000, 3, 3)

# sigma and xi nonstationary
gev_models[[7]]$parnames <- c("mu", "sigma0", "sigma1", "xi0", "xi1")
gev_models[[7]]$bound_lower <- c(0, 0, -200, -3, -3)
gev_models[[7]]$bound_upper <- c(5000, 1000, 200, 3, 3)

# all nonstationary
gev_models[[8]]$parnames <- c("mu0", "mu1", "sigma0", "sigma1", "xi0", "xi1")
gev_models[[8]]$bound_lower <- c(0, -500, 0, -200, -3, -3)
gev_models[[8]]$bound_upper <- c(5000, 500, 1000, 200, 3, 3)
