# Estimate alpha, beta and theta given lambda
estimateAlphaBetaTheta <- function(lambda, data, model, alpha.beta, theta) {
  
  #
  # Initial stuff
  #
  
  # Extract data
  y <- data$Cases
  r <- 2 - as.numeric(data$Reported)
  b <- data$b
  g <- data$g
  
  # Extract matrices
  B <- model$matrices$B
  X <- model$matrices$X
  BX <- model$matrices$BX
  Pt <- model$matrices$Pt
  Pd <- model$matrices$Pd
  Du <- model$matrices$Du
  Ps <- model$matrices$Ps
  Pw <- model$matrices$Pw
  
  # Extract fixed smoothing parameters
  kappa.u <- model$kappa["kappa.u"]
  kappa.b <- model$kappa["kappa.b"]
  kappa.w <- model$kappa["kappa.w"]
  kappa.s <- model$kappa["kappa.s"]
  
  # Get number of coefficients
  Ks <- ncol(B)
  Kw <- ncol(X)
  
  #
  # Estimate parameters
  #
  
  # Initial log-likelihood for the negative binomial distribution
  ll <- 10; ll.old <- 1
  it <- 1
  # Penalized iterative weighted least squares algorithm
  while (abs((ll - ll.old)/ll.old) > 1e-2 & it <= 50) {
    # Get alpha and beta
    alpha <- alpha.beta[1:Ks]
    beta  <- alpha.beta[(Ks + 1):(Ks + Kw)]
    # Linear predictor
    eta.s <- as.vector(B %*% alpha)
    eta.w <- as.vector(X %*% beta)
    eta <- eta.s + eta.w
    # Inverse link function: exp(eta)
    mu <- exp(eta)
    # Update theta
    opt.theta <- optim(
      par = log(theta),
      fn = function(log.theta) -sum(r*dnbinom(x = y, mu = mu, size = exp(log.theta), log = TRUE)),
      method = "BFGS",
      control = list(reltol = 1e-2),
      hessian = TRUE)
    theta <- exp(opt.theta$par)
    # Weights: W = [1/Var(y)]*dmu.deta^2
    W <- (1/(mu + mu^2/theta))*mu^2
    # Working variable: z = eta + (y - mu)*(1/dmu.deta)
    z <- eta + (y - mu)*(1/mu)
    # Unimodal constraint
    vu <- as.numeric(Du %*% alpha >= 0)
    Vu <- Diagonal(x = vu)
    Pu <- t(Du) %*% Vu %*% Du
    # Boundary constraint
    vb <- b*(eta.s >= g)
    Vb <- Diagonal(x = vb)
    XVbX <- t(B) %*% Vb %*% B
    XVbg <- t(B) %*% Vb %*% g
    # Penalty matrix
    P <- bdiag(lambda[1]*Pt + lambda[2]*Pd + kappa.u*Pu + kappa.b*XVbX + kappa.s*Ps, kappa.w*Pw)
    # Normal equations for weighted least squares
    XWX <- t(BX) %*% Diagonal(x = W*r) %*% BX
    XWz <- t(BX) %*% (W*r*z)
    XWX.P.inv <- solve(XWX + P)
    # Update alpha and beta
    alpha.beta <- XWX.P.inv %*% (XWz + c(kappa.b*as.vector(XVbg), rep(0, Kw)))
    # Calculate log-likelihood for given theta
    ll.old <- ll
    ll <- sum(r*dnbinom(x = y, mu = mu, size = theta, log = TRUE))
    # Update iterator
    it <- it + 1
  }
  
  #
  # Calculate information criterion
  #
  
  # Calculate effective dimension
  ed <- sum(diag(XWX.P.inv %*% XWX))
  # Calculate BIC
  bic <- -2*ll + log(sum(r))*ed
  
  #
  # Return output
  #
  
  return(list(
    alpha.beta = alpha.beta, alpha.beta.cov = XWX.P.inv,
    theta = theta,
    bic = bic))
}
