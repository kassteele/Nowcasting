genPriorDelayDist <- function(mean.delay, max.delay, p = 0.999) {
  #
  # Generate a prior delay distribution
  #
  # Description
  # This function generates a PMF for the prior delay distribution
  #
  # Arguments
  # mean.delay  Assumed mean delay in days
  # max.delay   Assumed maximum delay in days, where a fraction p of all cases have been reported
  #             Also the number of days back for which to make a delay correction
  # p           Fraction of reported cases within max.delay days. Default 99.9 %
  #
  # Value
  # PMF as a vector of length max.delay + 1
  # 
  # Details
  # Prior delay distribution is assumed to be Negative Binomial
  # Note that this prior delay distribution is disease specific!
  
  # Find the overdispersion parameter based on mean.delay, max.delay and p
  theta.delay <- exp(uniroot(
    f = function(x) qnbinom(p = p, mu = mean.delay, size = exp(x)) - max.delay,
    interval = c(0, 5))$root)

  # We expect 1 case on day 1 of the outbreak
  # log(f.priordelay) is then the boundary constraint for the trend surface
  f.priordelay <- 1*dnbinom(x = 0:max.delay, mu = mean.delay, size = theta.delay)
  f.priordelay <- f.priordelay/sum(f.priordelay)
  
  # Return output
  return(f.priordelay)
}
