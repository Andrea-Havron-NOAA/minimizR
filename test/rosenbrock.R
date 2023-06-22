library(minimizR)


objective <- function(x) {
  return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
}


gradient <- function(x) {
  return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]), 200 * (x[2] - x[1] * x[1]) ) )
}


x<-c(-1.2, 0)

opt<-minimizR(x,objective, gradient, control = list(tolerance = 1e-8, lb = c(-1.5, -1.5), ub = c(1.5,1.5)))

print(opt)
