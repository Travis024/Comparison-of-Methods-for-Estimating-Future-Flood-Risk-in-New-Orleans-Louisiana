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

# TODO -- here now!

##=============================================================================
## End
##=============================================================================
