##=============================================================================
## comparison.R
##
##=============================================================================

## Read Tony's calibration results
load("../output/mcmc_output_processed_normalgamma_12Jun2019.RData")
parameters_normalgamma <- parameters_posterior
load("../output/mcmc_output_processed_uniform_13Jun2019.RData")
parameters_uniform <- parameters_posterior
nmodel <- length(parameters_posterior)

## Read Mingxuan's calibration results
parameters_check <- vector("list", nmodel)
nparameter <- rep(NA, nmodel)
for (mm in 1:nmodel) {
  nparameter[mm] <- ncol(parameters_posterior[[mm]])
  parameters_check[[mm]] <- read.table(paste("../parameters_data/Data\ Frames/params",mm,".csv",sep=""), header = TRUE, sep=',')
  parameters_check[[mm]] <- parameters_check[[mm]][,-1]
}

## Check consistent parameter names on indexing the 8 models
for (mm in 1:nmodel) {
  print("===============================")
  print(paste("model ",mm,sep=""))
  print(paste("TW: ", colnames(parameters_posterior[[mm]]), sep=""))
  print(paste("MZ: ", colnames(parameters_check[[mm]]), sep=""))
}

## Swap Mingxuan's models 5 and 7
parameters_save <- parameters_check[[7]]
parameters_check[[7]] <- parameters_check[[5]]
parameters_check[[5]] <- parameters_save

## Fit KDEs and plot comparisons
density_normalgamma <- vector("list", nmodel)
density_uniform <- vector("list", nmodel)
density_check <- vector("list", nmodel)

for (mm in 1:nmodel) {
  density_normalgamma[[mm]] <- vector("list", ncol(parameters_normalgamma[[mm]]))
  density_uniform[[mm]] <- vector("list", ncol(parameters_uniform[[mm]]))
  density_check[[mm]] <- vector("list", ncol(parameters_posterior[[mm]]))
  for (pp in 1:nparameter[mm]) {
    density_normalgamma[[mm]][[pp]] <- density(parameters_normalgamma[[mm]][,pp])
    density_uniform[[mm]][[pp]] <- density(parameters_uniform[[mm]][,pp])
    density_check[[mm]][[pp]] <- density(parameters_check[[mm]][,pp])
  }
}

## Make plots
comp_plot <- function(mm, pp) {
  xlims <- quantile(c(density_check[[mm]][[pp]]$x, 
                      density_normalgamma[[mm]][[pp]]$x, 
                      density_uniform[[mm]][[pp]]$x), c(0,1))
  xlims[2] <- xlims[2]+0.15*diff(xlims)
  ylims <- quantile(c(density_check[[mm]][[pp]]$y, 
                      density_normalgamma[[mm]][[pp]]$y, 
                      density_uniform[[mm]][[pp]]$y), c(0,1))
  plot(density_check[[mm]][[pp]]$x, density_check[[mm]][[pp]]$y, type='l', col="coral", lwd=2,
       xlim=xlims, ylim=ylims, xlab=colnames(parameters_posterior[[mm]])[pp], ylab="Density")
  lines(density_normalgamma[[mm]][[pp]]$x, density_normalgamma[[mm]][[pp]]$y, type='l', col="steelblue", lty=1, lwd=2)
  lines(density_uniform[[mm]][[pp]]$x, density_uniform[[mm]][[pp]]$y, type='l', col="steelblue", lty=2, lwd=2)
  legend(xlims[1]+0.45*diff(xlims), ylims[2], bty="n", cex=0.62,
         c("Mingxuan","Tony, normal/gamma", "Tony, uniform"), lty=c(1,1,2), lwd=c(2,2,2), col=c("coral","steelblue","steelblue"))
}

## Example:
comp_plot(1,1)

##=============================================================================
## Check the return levels in 2016 and 2065 using each set of priors
##=============================================================================

## read temperature forcing
source("read_temperature_data.R")

returnperiod_proj <- 100
year_proj <- 2065
temp_proj <- temperature_forc[which(time_forc==year_proj)]

# model 1 (stationary)
params <- parameters_normalgamma[[1]]
n_ensemble <- nrow(parameters_normalgamma[[1]])
x100_stat <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=params[ii,"mu"], scale=params[ii,"sigma"], shape=params[ii,"xi"])})

# model 2 (mu nonstationary)
params <- parameters_normalgamma[[2]]
n_ensemble <- nrow(parameters_normalgamma[[2]])
x100_NSmu <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=(params[ii,"mu0"]+params[ii,"mu1"]*temp_proj), scale=params[ii,"sigma"], shape=params[ii,"xi"])})

# model 3 (sigma nonstationary)
params <- parameters_normalgamma[[3]]
n_ensemble <- nrow(parameters_normalgamma[[3]])
x100_NSsig <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=params[ii,"mu"], scale=exp(params[ii,"sigma0"]+params[ii,"sigma1"]*temp_proj), shape=params[ii,"xi"])})

# model 6 (mu and xi nonstationary)
params <- parameters_normalgamma[[6]]
n_ensemble <- nrow(parameters_normalgamma[[6]])
x100_NSmuxi <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=(params[ii,"mu0"]+params[ii,"mu1"]*temp_proj), scale=params[ii,"sigma"], shape=(params[ii,"xi0"]+params[ii,"xi1"]*temp_proj))})

#
lb <- -1000
ub <- 2e4
kde_stat <- density(x100_stat, from=lb, to=ub)
kde_NSmu <- density(x100_NSmu, from=lb, to=ub)
kde_NSsig <- density(x100_NSsig, from=lb, to=ub)
kde_NSmuxi <- density(x100_NSmuxi, from=lb, to=ub)

plot(kde_stat$x, kde_stat$y, type='l', lwd=2, col="black", xlab="Surge height [mm]", ylab="Density")
lines(kde_NSmu$x, kde_NSmu$y, col="steelblue", lty=2, lwd=2)
lines(kde_NSsig$x, kde_NSsig$y, col="coral", lty=2, lwd=2)
lines(kde_NSmuxi$x, kde_NSmuxi$y, col="seagreen", lty=2, lwd=2)


##=============================================================================
## with end year changing

returnperiod_proj <- 100

params <- parameters_normalgamma[[1]]
n_ensemble <- nrow(parameters_normalgamma[[1]])
x100_stat <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=params[ii,"mu"], scale=params[ii,"sigma"], shape=params[ii,"xi"])})

# model 2 (mu nonstationary)
params <- parameters_normalgamma[[2]]
n_ensemble <- nrow(parameters_normalgamma[[2]])

temp_proj <- temperature_forc[which(time_forc==2010)]
x100_2010 <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=(params[ii,"mu0"]+params[ii,"mu1"]*temp_proj), scale=params[ii,"sigma"], shape=params[ii,"xi"])})
temp_proj <- temperature_forc[which(time_forc==2020)]
x100_2020 <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=(params[ii,"mu0"]+params[ii,"mu1"]*temp_proj), scale=params[ii,"sigma"], shape=params[ii,"xi"])})
temp_proj <- temperature_forc[which(time_forc==2030)]
x100_2030 <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=(params[ii,"mu0"]+params[ii,"mu1"]*temp_proj), scale=params[ii,"sigma"], shape=params[ii,"xi"])})
temp_proj <- temperature_forc[which(time_forc==2040)]
x100_2040 <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=(params[ii,"mu0"]+params[ii,"mu1"]*temp_proj), scale=params[ii,"sigma"], shape=params[ii,"xi"])})
temp_proj <- temperature_forc[which(time_forc==2050)]
x100_2050 <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=(params[ii,"mu0"]+params[ii,"mu1"]*temp_proj), scale=params[ii,"sigma"], shape=params[ii,"xi"])})
temp_proj <- temperature_forc[which(time_forc==2060)]
x100_2060 <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=(params[ii,"mu0"]+params[ii,"mu1"]*temp_proj), scale=params[ii,"sigma"], shape=params[ii,"xi"])})
temp_proj <- temperature_forc[which(time_forc==2070)]
x100_2070 <- sapply(1:n_ensemble, function(ii) {qevd(1-1/returnperiod_proj, loc=(params[ii,"mu0"]+params[ii,"mu1"]*temp_proj), scale=params[ii,"sigma"], shape=params[ii,"xi"])})

lb <- -1000
ub <- 2e4
kde_stat <- density(x100_stat, from=lb, to=ub)
kde_2010 <- density(x100_2010, from=lb, to=ub)
kde_2020 <- density(x100_2020, from=lb, to=ub)
kde_2030 <- density(x100_2030, from=lb, to=ub)
kde_2040 <- density(x100_2040, from=lb, to=ub)
kde_2050 <- density(x100_2050, from=lb, to=ub)
kde_2060 <- density(x100_2060, from=lb, to=ub)
kde_2070 <- density(x100_2070, from=lb, to=ub)

plot(kde_stat$x, kde_stat$y, type='l', lwd=2, col="black", xlab="Surge height [mm]", ylab="Density", xlim=c(0,8000))
lines(kde_2010$x, kde_2010$y, col="aquamarine2", lty=2, lwd=2)
lines(kde_2030$x, kde_2030$y, col="aquamarine4", lty=2, lwd=2)
lines(kde_2050$x, kde_2050$y, col="steelblue1", lty=2, lwd=2)
lines(kde_2060$x, kde_2070$y, col="steelblue3", lty=2, lwd=2)

quantile(x100_stat, c(.5,.05,.95))
quantile(x100_2010, c(.5,.05,.95))
quantile(x100_2030, c(.5,.05,.95))
quantile(x100_2050, c(.5,.05,.95))

# TODO -- here now!

##=============================================================================
## End
##=============================================================================
