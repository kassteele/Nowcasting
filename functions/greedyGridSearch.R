greedyGridSearch <- function(fn, lower, upper, n.grid, start, logscale, ...) {
  # 
  # Optimization over a parameter grid
  # 
  # Description
  # This function does a greedy grid search in any dimension
  #
  # Arguments
  # fn	      A function to be minimized, with first argument the vector of parameters
  #           over which minimization is to take place. It should return a scalar result
  # lower	    Numeric vector containing the lower bounds on the parameter grid
  # upper     Numeric vector containing the upper bounds on the parameter grid
  # n.grid    Integer number determining grid length in every dimension
  # start     Optional numeric vector containing initial values for the parameters to be optimized over
  # logscale  Logical vector. If TRUE, a logarithmic scale is used for that parameter. It defaults to FALSE, i.e., a linear scale
  # ...       Further arguments to be passed to fn
  #
  # Value
  # A list with components:
  # par	      The best set of parameters found
  # value	    The value of fn corresponding to par
  # counts    A integer giving the number of calls to fn
  #
  # Details
  # A greedy algorithm is an algorithmic paradigm that follows the problem solving heuristic of making the
  # locally optimal choice at each stage with the hope of finding a global optimum. In many problems,
  # a greedy strategy works well if there are no local optima.
  
  # Get dimension of parameter space
  n.par <- length(lower)
  # Apply log-transformation to elements of lower and upper?
  lower <- ifelse(test = logscale, yes = log10(lower), no = lower)
  upper <- ifelse(test = logscale, yes = log10(upper), no = upper)
  # Set all possibles values on parameter grid. Result is n.grid x n.par matrix
  par.grid <- mapply(FUN = seq, from = lower, to = upper, length = n.grid)
  # Get initial grid index value
  if (missing(start)) {
    # If start is not given, take index half way the grid
    index <- rep(floor(n.grid/2), times = n.par)
  } else {
    # Else, take index of par.grid value closest to start
    start <- ifelse(test = logscale, yes = log10(start), no = start)
    index <- apply(X = abs(t(par.grid) - start), MARGIN = 1, FUN = which.min)
  }
  # Apply backtransformation to vectors of par.grid?
  par.grid <- matrix(
    mapply(FUN = function(x, logscale) ifelse(test = logscale, yes = 10^x, no = x), t(par.grid), logscale),
    ncol = n.par, byrow = TRUE)
  # Set initial parameter values
  par <- diag(par.grid[index, ])
  # Create n.par dimensional array of size rep(n.grid, n.par) filled with NA
  # Will be filled with function evaluations
  f.eval <- array(NA, dim = rep(n.grid, times = n.par))
  # Initial function value
  f.min <- Inf
  # Initital number of function evaluations
  n.eval <- 0
  # Set move to TRUE to enable search
  move <- TRUE
  while (move) {
    # Stop searching when no improvement
    move <- FALSE
    # For each parameter
    for (i in 1:n.par) {
      # Set candidate index vector to current index vector
      index.can <- index
      # Move i-th index one down and up
      for (j in max(index[i] - 1, 1):min(index[i] + 1, n.grid)) {
        # Set i-th index of index.can to j
        index.can[i] <- j
        # If there is no function value yet
        if (is.na(f.eval[t(index.can)])) {
          # Copy current parameters to candidates
          par.can <- par
          # Replace i-th candidate parameter by par.grid[j, i]
          par.can[i] <- par.grid[j, i]
          # Evaluate function for par.can
          f.new <- fn(par.can, ...)
          # Increase number of function evaluations by one
          n.eval <- n.eval + 1
          # Replace NA in f.eval by f.new
          f.eval[t(index.can)] <- f.new
        } else {
          # If there was already a function value, set f.new to that value
          f.new <- f.eval[t(index.can)]
        }
        # If f.new is smaller than f.min
        if (f.new < f.min) {
          # Update f.min by f.new
          f.min <- f.new
          # Update index
          index[i] <- j
          # Update current parameter values
          par <- par.can
          # Set move = TRUE to continue searching
          move <- TRUE
        }
      }
    }
  }
  # Return output
  return(list(par = par, value = f.min, counts = n.eval))
}
