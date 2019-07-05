#===============================================================================
# compute_weights.R
#
# This file reads in the output files from the marginal likelihood estimator
# for each model and combines them into a combined .RDS file with bma weights
# and marginal likelihoods.
#
# Requires the RData files from the 'bridge_sample.R' routine as input
#
# Original code: Vivek Srikrishnan (Penn State) 2017
# Modified code: Tony Wong (CU Boulder) 2018
#===============================================================================

rm(list=ls())

library(Hmisc)

types.of.priors <- 'uniform'
filename.likelihood <- paste('log_marginal_likelihood_',types.of.priors,'.rds', sep='')
filename.weights <- paste('bma_weights_',types.of.priors,'.rds', sep='')

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/students/Travis/nola_surge/R')
  path.ml <- paste('/Users/tony/codes/students/Travis/nola_surge/output/bma/',types.of.priors,'/', sep="")
  path.out <- '/Users/tony/codes/students/Travis/nola_surge/output'
  path.R <- '/Users/tony/codes/students/Travis/nola_surge/R'
} else {
  # ???
}

source(paste(path.R,'parameter_setup.R',sep='/'))
source(paste(path.R,'read_temperature_data.R',sep='/'))

bma.weights <- rep(NA, nmodel)
log.marg.lik <- rep(NA, nmodel)

files <- list.files(path=path.ml, full.names=TRUE, recursive=FALSE)

for (file in files) {
  load(file)
  ##site <- capitalize(toString(station))
  gev.model <- as.numeric(substr(unlist(strsplit(file, split='_'))[4],1,1))
  log.marg.lik[gev.model] <- ml[length(ml)]
}

ml <- log.marg.lik
ml.scale <- ml - max(ml,na.rm=TRUE)
for (gev.model in 1:nmodel) {
  if (!is.na(log.marg.lik[gev.model])) {
    bma.weights[gev.model] <- exp(ml.scale[gev.model])/sum(exp(ml.scale), na.rm=TRUE)
  }
}

saveRDS(log.marg.lik, paste(path.out,filename.likelihood, sep="/"))
saveRDS(bma.weights, paste(path.out,filename.weights, sep="/"))

#===============================================================================
# End
#===============================================================================
