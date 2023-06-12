library(minimizR)


objective <- function(x) {
  return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
}


gradient <- function(x) {
  return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]), 200 * (x[2] - x[1] * x[1]) ) )
}


x<-c(-1.2, 0)

opt<-minimize(x,objective, gradient, control = list(tolerance = 1e-8, minb = c(-1.5, -1.5), maxb = c(1.5 ,1.5)))
#opt <- optim(x, objective, gradient, method = "BFGS", control= list(trace = 2))
print(opt)
