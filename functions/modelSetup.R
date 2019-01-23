modelSetup <- function(data, ord = 2, kappa = list(u = 1e6, b = 1e6, w = 0.01, s = 1e-6)) {
  #
  # Create model setup
  #
  # Description
  # Function that sets up the nowcasting model
  #
  # Arguments
  # data  Dataframe, output from dataSetup
  #
  # Value
  # List with:
  # matrices  List of model matrices and penalty matrices
  # kappa     Vector with fixed smoothing parameters for constraints
  
  #
  # Initial stuff
  #
  
  # Filter on records with Est == 1
  data <- data %>% filter(Est == 1)
  
  # Extract dimensions
  t <- unique(data$t)
  d <- unique(data$d)
  T  <- length(t)
  D1 <- length(d)
  
  #
  # Model matrices
  #
  
  # B-spline basis matrix for smooth surface
  Bt <- bbase(x = t, k = max(4, floor( T/5)), deg = 3)
  Bd <- bbase(x = d, k = max(4, floor(D1/5)), deg = 3)
  B <- kronecker(Bd, Bt)
  
  # Model matrix for weekday effect
  # Because intercept is included in B-spline basis, drop first column (Monday = reference) of X
  X <- sparse.model.matrix(~ Day, data = data)[, -1]
  
  # cbind them together
  BX <- cbind(B, X)
  
  # Get number of coefficients
  Kt <- ncol(Bt)
  Kd <- ncol(Bd)
  Kw <- ncol(X)
  
  #
  # Difference operator and penalty matrices
  #
  
  # Difference operator matrices
  Dt <- kronecker(Diagonal(Kd), diff(Diagonal(Kt), diff = ord)) # Smoothness in t direction
  Dd <- kronecker(diff(Diagonal(Kd), diff = 2), Diagonal(Kt))   # Smoothness in d direction
  Du <- kronecker(diff(Diagonal(Kd), diff = 2), Diagonal(Kt))   # Unimodal in d direction
  
  # Penalty matrices
  Pt <- t(Dt) %*% Dt
  Pd <- t(Dd) %*% Dd
  Pw <- Diagonal(Kw)
  Ps <- Diagonal(Kt*Kd)
  
  #
  # Fixed smoothing parameters
  #
  
  # kappa.u and kappa.b are large for asymmetric penalty
  # kappa.s is very small for ridge penalty on surface
  # kappa.w is small for ridge penalty on weekday effect
  kappa.u <- kappa$u
  kappa.b <- kappa$b
  kappa.w <- kappa$w
  kappa.s <- kappa$s
  
  #
  # Return output
  #

  return(list(
    matrices = list(B = B, X = X, BX = BX, Pt = Pt, Pd = Pd, Du = Du, Ps = Ps, Pw = Pw),
    kappa = c(kappa.u = kappa.u, kappa.b = kappa.b, kappa.w = kappa.w, kappa.s = kappa.s)))
  
}
