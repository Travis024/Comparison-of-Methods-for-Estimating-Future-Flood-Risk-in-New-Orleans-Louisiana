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
.NP.deoptim <- 100      # number of DE population members (at least 10*[# parameters])
.niter.deoptim <- 100   # number of DE iterations
output.dir <- '../output/'
data.dir <- '../parameters_data/'
DO_PROCESSING <- TRUE  # TRUE if you need to run the processing
  # TRUE -> detrend with moving window 
  # FALSE -> read in some previous RDS processing results,
  # with the filenames defined below:
filename.many <- ''
filename.nola <- ''
threshold_missing_data <- 0.9 # filter out years that are missing more data than this (%)
fillvalue <- -32767

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
## helper functions...

# for processing data
source("process_data_function.R")

# log-likelihood function, prior and posterior
source("likelihood_function.R")

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

TODO
##=============================================================================


##=============================================================================
## Define likelihood, priors and posterior
## priors -- based on long tide gauges as in Wong et al 2018 (doi: ...)



##=============================================================================
## Calibration



##=============================================================================
## Burn-in and thinning



##=============================================================================
## End
##=============================================================================
