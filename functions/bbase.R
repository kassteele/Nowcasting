bbase <- function(x, x.min = min(x), x.max = max(x), k = 15, deg = 3, sparse = TRUE) {
  #
  # bbase
  #
  # Description
  # Generates design matrix for B-splines
  #
  # Arguments
  # x       A numeric vector of values at which to evaluate the B-spline functions 
  # x.min   Lowest value, min(x)
  # x.max   Highest value, max(x)
  # k       Number of B-spline basis functions
  # deg     Degree of the B-spline basis function. Default is cubic B-splines
  # sparse  Logical indicating if the result should inherit from class "sparseMatrix" (from package Matrix)
  #
  # Value
  # Matrix B-spline basis functions
  
  dx <- (x.max - x.min)/(k - deg)
  knots <- seq(from = x.min - deg*dx, to = x.max + deg*dx, by = dx)
  B <- splines::splineDesign(x = x, knots = knots, ord = deg + 1, outer.ok = TRUE, sparse = sparse)
  return(B)
}
