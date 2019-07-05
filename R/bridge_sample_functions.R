#===============================================================================
# This file contains the working script for estimating the          #
# marginal likelihood of different GPD surge models and data        #
# experiments.                                                      #
#                                                                   #
# This revised version calculates all likelihood estimates as a parallel
# 'foreach' loop. (TODO)
#===============================================================================


#===============================================================================
# compute the bridge sampling estimate using the iterative procedure from Meng and Wong (1996)
# starting value is not quite the reciprocal importance sampling estimate
# from Gelfand and Dey (1994)
# due to infinite values obtained when exponentiating, we start with the reciprocal mean of
recip.imp.samp <- function(log.p,log.imp) {
  log.ratio <- log.imp - log.p
  -log(mean(exp(log.ratio)))
}
#===============================================================================


#===============================================================================
# function to update the bridge sampling estimate at each iteration
# norm.const is the log normalizing constant estimate from the previous iteration
# post, imp are lists with log.p and log.imp passing the associated log-likelihoods
bridge.samp.iter <- function(log.norm.const,
                        post,
                        imp) {

  # normalize posterior likelihoods based on previous normalizing constant estimate
  # some samples (mostly importance samples) might have infinite posterior log-likelihoods
  post.log.p.norm <- post$log.p[is.finite(post$log.p)] - log.norm.const
  imp.log.p.norm <- imp$log.p[is.finite(imp$log.p)] - log.norm.const

  post.log.imp <- post$log.imp[is.finite(post$log.p)]
  imp.log.imp <- imp$log.imp[is.finite(imp$log.p)]

  # get number of samples
  post.num <- length(post.log.p.norm)
  imp.num <- length(imp.log.p.norm)

  # compute updated estimate numerator and denominator
  imp.mean <- mean(exp(imp.log.p.norm)/(imp.num*exp(imp.log.imp)+post.num*exp(imp.log.p.norm)))
  post.mean <- mean(exp(post.log.imp)/(imp.num*exp(post.log.imp)+post.num*exp(post.log.p.norm)))

  # return updated estimate
  log.norm.const + log(imp.mean) - log(post.mean)
}
#===============================================================================


#===============================================================================
# compute the relative standard error of the bridge sampling estimator
# we can treat the posterior samples as iid due to re-sampling from the posterior,
# so we use the error formula from Fruhwirth-Schnatter (2004) with the spectral density
# at frequency 0 set equal to 1.
bridge.samp.rel.err <- function(log.norm.const,
                                post,
                                imp) {

  # normalize posterior likelihoods based on previous normalizing constant estimate
  # some samples (mostly importance samples) might have infinite posterior log-likelihoods
  post.log.p.norm <- post$log.p[is.finite(post$log.p)] - log.norm.const
  imp.log.p.norm <- imp$log.p[is.finite(imp$log.p)] - log.norm.const

  post.log.imp <- post$log.imp[is.finite(post$log.p)]
  imp.log.imp <- imp$log.imp[is.finite(imp$log.p)]

  # get number of samples
  post.num <- length(post.log.p.norm)
  imp.num <- length(imp.log.p.norm)
  imp.factor <- imp.num/(imp.num+post.num)
  post.factor <- post.num/(imp.num+post.num)

  # compute bridging function process estimates for both sequences
  post.bridge.int <- exp(post.log.imp)/(imp.factor*exp(post.log.imp)+post.factor*(exp(post.log.p.norm)))
  imp.bridge.int <- exp(imp.log.p.norm)/(imp.factor*exp(imp.log.imp)+post.factor*(exp(imp.log.p.norm)))

  # return squared relative error estimate
  (var(imp.bridge.int)/mean(imp.bridge.int)^2)/imp.num + (var(post.bridge.int)/mean(post.bridge.int)^2)/post.num
}
#===============================================================================
