#===============================================================================
# This file contains the working script for estimating the marginal likelihood
# of different GPD surge models and data experiments.
#
# This revised version calculates all likelihood estimates as a parallel
# 'foreach' loop.
#
# Original code: Vivek Srikrishnan (Penn State) 2017
# Modified code: Tony Wong (CU Boulder) 2018
#===============================================================================

# import libraries
library(mvtnorm)  # used for multivariate normal samples and densities
#library(coda)     # used to compute spectral density of MCMC draws at frequency 0
library(extRemes) # used to compute GEV densities within posterior likelihood function
library(foreach)
library(doParallel)
library(ncdf4)

calib_date <- '13Jun2019'
type.of.priors <- 'uniform'     # can be 'uniform' or 'normalgamma'

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  path.R <- '/Users/tony/codes/students/Travis/nola_surge/R'
  setwd(path.R)
  # set data and save directories
  path.data <- '/Users/tony/codes/students/Travis/nola_surge/parameters_data'
  path.output <- '/Users/tony/codes/students/Travis/nola_surge/output'
  path.save <- paste('/Users/tony/codes/students/Travis/nola_surge/output/bma/',type.of.priors,'/', sep='')
  nnode <- 1          # number of CPUs to use
} else {
  # ???
}

# import file containing the log likelihood calculations
source(paste(path.R,'likelihood.R',sep='/'))

# read temperature covariate data
source(paste(path.R,'read_temperature_data.R',sep='/'))

# set up parameters
source(paste(path.R,'parameter_setup.R',sep='/'))

# list to store the actual output
output <- vector('list', nmodel)  # named just by number
gev.models <- seq(from=1, to=nmodel)

# data frame to store the experimental details
experiments <- expand.grid(station  ="grandisle",
                           gev.model=gev.models)
n_experiments <- nrow(experiments)
output <- vector('list', n_experiments)

# data from Grand Isle
threshold_missing_data <- 0.9 # filter out years that are missing more data than this (%)
fillvalue <- -32767
source(paste(path.R,'data_processing.R',sep='/'))
file.nola <- paste(path.data,"h765a_grandisle.csv", sep="/")
data.nola <- process_tg_data(file.nola)

#cores = detectCores()
#cl <- makeCluster(cores[1]-1) #not to overload your computer
cl <- makeCluster(nnode)
print(paste('Starting cluster with ',nnode,' cores', sep=''))
registerDoParallel(cl)

source('bridge_sample_functions.R')

export.names <- c('bridge.samp.rel.err','bridge.samp.iter','recip.imp.samp','experiments','log_post_gev','log_like_gev','log_prior_gev','path.R','calib_date','type.of.priors')

finalOutput <- foreach(ee=1:n_experiments,
                            .packages=c('mvtnorm','extRemes','ncdf4'),
                            .export=export.names,
                            .inorder=FALSE) %dopar% {

  setwd(path.R)
  source(paste(path.R,'likelihood.R',sep='/'))
  source(paste(path.R,'read_temperature_data.R',sep='/'))
  # get parameters for this particular experiment
  print(experiments[ee,])
  station <- experiments[ee,'station']
  gev.model <- experiments[ee,'gev.model']

  # set output (saved as .RData; to be collected into a single output file later) file name
  filename.out <- paste('ml_',station,'_',gev.model,'.RData',sep='')
  if (file.exists(paste(path.save, filename.out, sep='/'))) {
     #stop('Output file already exists!')
     print('Output file already exists!')
     output[[ee]] <- 'done!'
  } else {

    # read in calibration output file
    print('loading calibration file...')

    setwd(path.output)
    filename.priors <- Sys.glob(paste('surge_priors_',type.of.priors,'_*','.RData',sep='')) # is in the output directory

    # use this if multiple files exist for the same location and prior
    load(filename.priors)
    if (exists('calib_date')) {
      setwd(path.output)
      filename.calib <- paste('mcmc_output_processed_',type.of.priors,'_',calib_date,'.RData',sep='')
    } else {
      setwd(path.output)
      filename.calib <- Sys.glob(paste('mcmc_output_processed_',type.of.priors,'_*','.RData',sep=''))
    }
    load(filename.calib)  # gives parameters_posterior[[m]]  m in 1:nmodel

    print('done!')

    # these chains are burned in and thinned, so use the whole thing.
    nsamples <- nrow(parameters_posterior[[gev.model]])

    # set number of samples to use for estimate
    post.samp.num <- nsamples
    imp.samp.num <- nsamples

    # burn in samples and log.p values
    post.samples <- parameters_posterior[[gev.model]]
    if (gev.model > 1) {aux <- trimmed_forcing(data.nola[,"year"], time_forc, temperature_forc)$temperature
    } else {aux <- NULL}
    post.ll <- apply(post.samples, 1, log_post_gev,
                            parnames=gev_models[[gev.model]]$parnames,
                            data_calib=data.nola[,"lsl_max"],
                            priors=priors[[gev.model]],
                            auxiliary=aux)

    # fit normal approximation to the posteriors
    post.mean <- colMeans(post.samples)
    post.cov <- cov(post.samples)

    # get posterior samples
    print('sampling from posterior distribution...')

    samp.names <- c('samples','log.imp','log.p')
    post.samp <- setNames(vector("list",length(samp.names)),samp.names)
    samp.idx <- sample(x=nrow(post.samples), size=post.samp.num, replace=TRUE)
    post.samp$samples <- post.samples
    post.samp$samples <- post.samples[samp.idx,]
    # get posterior log-likelihood of sampled posterior values
    post.samp$log.p <- post.ll
    post.samp$log.p <- post.ll[samp.idx]
    # get importance log-likelhood of posterior samples
    post.samp$log.imp <- dmvnorm(x=post.samp$samples, mean=post.mean, sigma=post.cov, log=TRUE)

    print('done!')

    # get importance samples and likelihood
    print('sampling from importance distribution...')

    imp.samp <- setNames(vector("list",length(samp.names)),samp.names)
    imp.samp$samples <- rmvnorm(n=imp.samp.num, mean=post.mean, sigma=post.cov)
    imp.samp$log.imp <- dmvnorm(x=imp.samp$samples, mean=post.mean, sigma=post.cov, log=TRUE)
    colnames(imp.samp$samples) <- colnames(post.samp$samples)
    # compute posterior log-likelihood of importance samples
    imp.samp$log.p <- apply(imp.samp$samples, 1, log_post_gev,
                            parnames=gev_models[[gev.model]]$parnames,
                            data_calib=data.nola[,"lsl_max"],
                            priors=priors[[gev.model]],
                            auxiliary=aux)

    print('done!')

    print('beginning bridge sampling recursion...')

    # set tolerance for halting of iteration
    TOL <- 1e-10

    # initialize storage for estimates
    ml <- mat.or.vec(nr=1,nc=1)

    # initialize with starting value
    # we can't quite start with the reciprocal importance sampling estimate from
    # Gelfand and Dey (1994) due to numerics (we get 0 values when we exponentiate
    # the difference of the importance log-densities and posterior log-likelihoods), so we just
    # average the ratios on a log scale.
    ml[1] <- -mean(post.samp$log.imp - post.samp$log.p)
    ml[2] <- bridge.samp.iter(ml[1], post.samp[c('log.p','log.imp')], imp.samp[c('log.p','log.imp')])

    # iterate until within tolerance.
    t <- 2
    while (abs(ml[t] - ml[t-1]) >= TOL) {
      ml[t+1] <- bridge.samp.iter(ml[t], post.samp[c('log.p', 'log.imp')], imp.samp[c('log.p', 'log.imp')])
      t <- t+1
    }

    print('done!')

    print('computing relative standard error of estimate')

    # compute the relative standard error of the bridge sampling estimator
    # we can treat the posterior samples as iid due to re-sampling from the posterior,
    # so we use the error formula from Fruhwirth-Schnatter (2004) with the spectral density
    # at frequency 0 set equal to 1.

    re.sq <- bridge.samp.rel.err(ml[length(ml)], post.samp[c('log.p','log.imp')], imp.samp[c('log.p','log.imp')])

    # save result of run
    # if save directory doesn't exist, create it
    #ifelse(!dir.exists(path.save), dir.create(path.save), FALSE)
    setwd(path.save)

    save(list=c('post.samp','imp.samp', 'ml', 're.sq', 'gev.model'), file=filename.out)
    output[[ee]] <- 'done!'
  }
}
stopCluster(cl)

#data_many <- finalOutput
#names(data_many) <- names(data_set)


#===============================================================================
# end
#===============================================================================
