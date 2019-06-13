##=============================================================================
## Calibration of (non-)stationary GEV models based on Grand Isle, LA tide
## gauge data.
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##=============================================================================



##=============================================================================
## Preliminary set-up
rm(list=ls())

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/students/Travis/Comparison-of-Methods-for-Estimating-Future-Flood-Risk-in-New-Orleans-Louisiana/R')
}
output.dir <- '../output/'
data.dir <- '../parameters_data/'

# Processing of data
DO_PROCESSING <- TRUE  # TRUE if you need to run the processing
  # TRUE -> detrend with moving window 
  # FALSE -> read in some previous RDS processing results,
  # with the filenames defined below:
filename.many <- ''
filename.nola <- ''
threshold_missing_data <- 0.9 # filter out years that are missing more data than this (%)
fillvalue <- -32767

# MLE optimization for priors
type.of.priors <- "uniform" # "normalgamma" or "uniform"
NP.deoptim <- 100      # number of DE population members (at least 10*[# parameters])
niter.deoptim <- 100   # number of DE iterations
F.deoptim <- 0.8
CR.deoptim <- 0.9

# Convergence and autocorrelation check
gr.max <- 1.1
lag.max <- 2000
acf.max <- 0.05

## Install/load relevant libraries
library(extRemes)
library(zoo)
library(date)
library(Hmisc)
library(adaptMCMC)
library(DEoptim)
library(ncdf4)
##=============================================================================


  
##=============================================================================
## helper functions and some set up

# for processing data
source("data_processing.R")

# log-likelihood function, prior and posterior
source("likelihood.R")

# set up parameters for all 8 candidate model structures
source("parameter_setup.R")
##=============================================================================



##=============================================================================
## Read Grand Isle tide gauge data, from UHSLC:
##    http://uhslc.soest.hawaii.edu/data/csv/rqds/atlantic/hourly/h765a.csv

file.nola <- paste(data.dir,"h765a_grandisle.csv", sep="")
data.nola <- process_tg_data(file.nola)
##=============================================================================



##=============================================================================
## Read other long data sets, also from UHSLC, to fit priors

data.dir.priors <- paste(data.dir, "tide_gauge_long/", sep="")
files.priors <- list.files(path=data.dir.priors, pattern="csv")
data.priors <- vector('list', length(files.priors))

for (dd in 1:length(files.priors)) {
  filename.new <- paste(data.dir.priors, files.priors[dd], sep='')
  data.new <- process_tg_data(filename.new)
  names(data.priors)[dd] <- substr(files.priors[dd], start=1, stop=7)
  data.priors[[dd]] <- data.new
}
##=============================================================================



##=============================================================================
## Read historical and future temperature data

source("read_temperature_data.R")
##=============================================================================



##=============================================================================
## Fit maximum likelihood estimates for the long tide gauge data sets, and
## fit prior distributions to these for all 8 model structures

# will want something like
# auxiliary <- trimmed_forcing(data.priors[[1]][,1], time_forc, temperature_forc)$temperature
# log_like_gev(parameters=c(200,1,100,1), parnames=gev_models[[2]]$parnames, data_calib=data.priors[[1]][,"lsl_max"], auxiliary=auxiliary)

deoptim.priors <- vector('list', nmodel)
for (mm in 1:nmodel) {
  deoptim.priors[[mm]] <- mat.or.vec(length(data.priors), length(gev_models[[mm]]$parnames))
  rownames(deoptim.priors[[mm]]) <- names(data.priors)
  colnames(deoptim.priors[[mm]]) <- gev_models[[mm]]$parnames
}

for (dd in 1:length(data.priors)) {
  print(paste('starting to calculate MLE GEV parameters for tide gauge data set ',dd,' / ',length(data.priors),sep=''))
  tbeg0 <- proc.time()
  # GEV model fitting
  for (mm in 1:nmodel) {
    print(paste('  - starting DE optimization for model ',mm,'...', sep=''))
    tbeg <- proc.time()
    # if tide gauge record starts before auxiliary forcing, clip it
    if(data.priors[[dd]][1,"year"] < time_forc[1]) {
      irem <- which(data.priors[[dd]][,"year"] < time_forc[1])
      data.priors[[dd]] <- data.priors[[dd]][-irem,]
    }
    # if tide gauge record ends after auxiliary forcing, clip it
    if(max(data.priors[[dd]][,"year"]) > max(time_forc)) {
      irem <- which(data.priors[[dd]][,"year"] > max(time_forc))
      data.priors[[dd]] <- data.priors[[dd]][-irem,]
    }
    # set up storm surge covariate
    if (mm > 1) {auxiliary <- trimmed_forcing(data.priors[[dd]][,"year"], time_forc, temperature_forc)$temperature
    } else {auxiliary <- NULL}
    # calibration
    out.deoptim <- DEoptim(neg_log_like_gev, lower=gev_models[[mm]]$bound_lower, upper=gev_models[[mm]]$bound_upper,
                           DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                           parnames=gev_models[[mm]]$parnames, data_calib=data.priors[[dd]][,"lsl_max"], auxiliary=auxiliary)
    deoptim.priors[[mm]][dd,] <- out.deoptim$optim$bestmem
    colnames(deoptim.priors[[mm]]) <- gev_models[[mm]]$parnames
    tend <- proc.time()
    print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  }
  tend0 <- proc.time()
  print(paste('... done. Took ',round(as.numeric(tend0-tbeg0)[3]/60,2),' minutes', sep=''))
}

## Fit priors -- based on long tide gauges as in Wong et al 2018 (doi: ...)

print('fitting prior distributions to the MLE parameters...')

# fit gamma and normal priors
# -> centered at the medians
# -> with standard deviation equal to half the max-min range
#    (or do empirical sd? might underestimate though - take wider)

# assign which parameters have which priors
if (type.of.priors=="normalgamma") {
  gamma.priors <- c('mu','mu0','sigma','sigma0')
  normal.priors <- c('mu1','sigma1','xi','xi0','xi1')
  uniform.priors <- NULL
} else if (type.of.priors=="uniform") {
  gamma.priors <- NULL
  normal.priors <- NULL
  uniform.priors <- c('mu','mu0','mu1','sigma','sigma0','sigma1','xi','xi0','xi1')
} else {print("ERROR: unknown type.of.priors")}

priors <- vector('list', nmodel)
for (mm in 1:nmodel) {
  priors[[mm]] <- vector('list', length(gev_models[[mm]]$parnames)); names(priors[[mm]]) <- gev_models[[mm]]$parnames
  for (par in gev_models[[mm]]$parnames) {
    priors[[mm]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
    if(!is.na(match(par, uniform.priors))) {
      names(priors[[mm]][[par]]) <- c('type','lower','upper'); priors[[mm]][[par]]$type <- 'uniform'
      priors[[mm]][[par]]$lower <- gev_models[[mm]]$bound_lower[match(par,gev_models[[mm]]$parnames)]
      priors[[mm]][[par]]$upper <- gev_models[[mm]]$bound_upper[match(par,gev_models[[mm]]$parnames)]
    } else if(!is.na(match(par, gamma.priors))) { # shape=alpha, rate=beta, mean=shape/rate, var=shape/rate^2
      names(priors[[mm]][[par]]) <- c('type','shape','rate'); priors[[mm]][[par]]$type <- 'gamma'
      priors[[mm]][[par]]$rate <- median(deoptim.priors[[mm]][,par]) / (0.5*(max(deoptim.priors[[mm]][,par])-min(deoptim.all[[mm]][,par])))^2
      priors[[mm]][[par]]$shape <- median(deoptim.priors[[mm]][,par]) * priors[[mm]][[par]]$rate
    } else if(!is.na(match(par, normal.priors))) {
      names(priors[[mm]][[par]]) <- c('type','mean','sd'); priors[[mm]][[par]]$type <- 'normal'
      priors[[mm]][[par]]$mean <- median(deoptim.priors[[mm]][,par])
      priors[[mm]][[par]]$sd   <- 0.5*(max(deoptim.priors[[mm]][,par])-min(deoptim.priors[[mm]][,par]))
      #priors[[mm]][[par]]$sd   <- sd(deoptim.all[[model]][,par])
    }
  }
}

print('...done.')
print(paste('saving priors and DE optim output as .rds files to read and use later...',sep=''))

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.priors <- paste(output.dir,'surge_priors_',type.of.priors,"_",today,'.RData', sep='')
save(list=c("priors","deoptim.priors"), file=filename.priors)
##=============================================================================



##=============================================================================
## Calibration for NOLA data

niter_mcmc <- 5e5 # note that this is enough for stationary model, but maybe not quite enough for 5- or 6-parameter nonstationary
nnode_mcmc <- 2   # with 2 chains, 1e5 iterations requires 106 sec on laptop
gamma_mcmc <- 0.66
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
amcmc_out <- vector("list", nmodel)
filename.mcmc <- paste(output.dir,'mcmc_output_',type.of.priors,"_",today,'.RData', sep='')

for (mm in 1:nmodel) {
  if (mm > 1) {auxiliary <- trimmed_forcing(data.nola[,"year"], time_forc, temperature_forc)$temperature
  } else {auxiliary <- NULL}
  accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(gev_models[[mm]]$parnames)
  step_mcmc <- 0.01*(gev_models[[mm]]$bound_upper-gev_models[[mm]]$bound_lower)
  
  # NOLA initial parameter estimates
  out.deoptim <- DEoptim(neg_log_like_gev, lower=gev_models[[mm]]$bound_lower, upper=gev_models[[mm]]$bound_upper,
                         DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                         parnames=gev_models[[mm]]$parnames, data_calib=data.nola[,"lsl_max"], auxiliary=auxiliary)
  initial_parameters <- out.deoptim$optim$bestmem
  
  # calibration by MCMC
  tbeg <- proc.time()
  if (nnode_mcmc==1) {
    amcmc_out[[mm]] <- MCMC(log_post_gev, n=niter_mcmc, init=initial_parameters, adapt=TRUE, acc.rate=accept_mcmc,
                            scale=step_mcmc, gamma=gamma_mcmc, list=TRUE, n.start=5000,
                            parnames=gev_models[[mm]]$parnames, data_calib=data.nola[,"lsl_max"],
                            priors=priors[[mm]], auxiliary=auxiliary)
  } else if (nnode_mcmc > 1) {
    amcmc_out[[mm]] <- MCMC.parallel(log_post_gev, n=niter_mcmc, init=initial_parameters,
                                     n.chain=nnode_mcmc, n.cpu=nnode_mcmc, packages='extRemes',
                                     adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc, gamma=gamma_mcmc, 
                                     list=TRUE, n.start=5000,
                                     parnames=gev_models[[mm]]$parnames, data_calib=data.nola[,"lsl_max"],
                                     priors=priors[[mm]], auxiliary=auxiliary)
  }
  tend <- proc.time()
  timer <- round((tend-tbeg)[3]/60,2)
  print(paste(nnode_mcmc," chains x ",niter_mcmc," iterations took ",timer," minutes",sep=""))
  print(paste('saving MCMC output as .RData file...',sep=''))
  save("amcmc_out", file=filename.mcmc)
}
##=============================================================================



##=============================================================================
## Burn-in and thinning

##=========##
## Burn-in ##
##=========##

niter.test <- seq(from=1e5, to=(niter_mcmc-1), by=1e5)
gr.test <- matrix(NA, nrow=nmodel, ncol=length(niter.test))

if(nnode_mcmc == 1) {
  # don't do GR stats, just cut off first half of chains
  print('only one chain; will lop off first half for burn-in instead of doing GR diagnostics')
} else if(nnode_mcmc > 1) {
  # this case is FAR more fun
  for (mm in 1:nmodel) {
    # accumulate the names of the soon-to-be mcmc objects
    string.mcmc.list <- 'mcmc1'
    for (cc in 2:nnode_mcmc) {
      string.mcmc.list <- paste(string.mcmc.list, ', mcmc', cc, sep='')
    }
    for (i in 1:length(niter.test)) {
      for (cc in 1:nnode_mcmc) {
        # convert each of the chains into mcmc object
        eval(parse(text=paste('mcmc',cc,' <- as.mcmc(amcmc_out[[mm]][[cc]]$samples[(niter.test[i]+1):niter_mcmc,])', sep='')))
      }
      eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))
      gr.test[mm,i] <- as.numeric(gelman.diag(mcmc_chain_list, autoburnin=FALSE)[2])
    }
  }
} else {print('error - n_node000 < 1 makes no sense')}

# hack off first ifirst iterations for burn in
ifirst <- NA
if (nnode_mcmc==1) {
  ifirst <- round(0.5*niter_mcmc)
} else {
  i <- 1
  CONVERGED <- all(gr.test[,i] < gr.max)
  while (!CONVERGED) {
    i <- i+1
    CONVERGED <- all(gr.test[,i] < gr.max)
  }
  ifirst <- niter.test[i]
}

chains_burned <- vector("list", nmodel)
if (nnode_mcmc > 1) {
  for (mm in 1:nmodel) {
    chains_burned[[mm]] <- vector('list', nnode_mcmc)
    for (cc in 1:nnode_mcmc) {
      chains_burned[[mm]][[cc]] <- amcmc_out[[mm]][[cc]]$samples[(ifirst+1):niter_mcmc,]
    }
  }
} else {
  for (mm in 1:nmodel) {
    chains_burned[[mm]] <- amcmc_out[[mm]]$samples[(ifirst+1):niter_mcmc,]
  }
}

##==========##
## Thinning ##
##==========##

maxlag <- 0

for (mm in 1:nmodel) {
  for (cc in 1:nnode_mcmc) {
    for (pp in 1:length(gev_models[[mm]]$parnames)) {
      acf_tmp <- acf(chains_burned[[mm]][[cc]][,pp], lag.max=lag.max, plot=FALSE)
      new <- acf_tmp$lag[which(acf_tmp$acf < acf.max)[1]]
      if (maxlag < new) {
        print(paste(mm,cc,pp,"Updating maxlag to",new))
        maxlag <- new
      }
    }
  }
}

chains_burned_thinned <- chains_burned # initialize
parameters_posterior <- vector("list", nmodel)
if (nnode_mcmc > 1) {
  for (mm in 1:nmodel) {
    for (cc in 1:nnode_mcmc) {
      chains_burned_thinned[[mm]][[cc]] <- chains_burned[[mm]][[cc]][seq(from=1, to=(niter_mcmc - ifirst), by=maxlag),]
    }
    parameters_posterior[[mm]] <- chains_burned_thinned[[mm]][[1]]
    for (cc in 2:nnode_mcmc) {
      parameters_posterior[[mm]] <- rbind(parameters_posterior[[mm]], chains_burned_thinned[[mm]][[cc]])
    }
  }
} else {
  for (mm in 1:nmodel) {
    chains_burned_thinned[[mm]] <- chains_burned[[mm]][seq(from=1, to=nrow(chains_burned), by=maxlag),]
    parameters_posterior[[mm]] <- chains_burned_thinned[[mm]]
  }
}

# Name the columns of the posterior parameters for reading easier later
for (mm in 1:nmodel) {
  colnames(parameters_posterior[[mm]]) <- gev_models[[mm]]$parnames
}
##=============================================================================



##=============================================================================
## Write parameters file

filename.posterior <- paste(output.dir,"mcmc_output_processed_",type.of.priors,"_",today,".RData", sep="")
save(parameters_posterior, file=filename.posterior)
##=============================================================================



##=============================================================================
## End
##=============================================================================
