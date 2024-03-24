
# minimizR
minimizR is a l-bfgs function minimizer that works with freely estimable or box constrained parameter sets.

# Install from R

```r
remotes::install_github("nmfs-fish-tools/minimizR")
library(minimizR)
```



## Freely Estimated Example

Here is an example of minimizing the Rosenbrock function with freely estimated parameters.

```r
library(minimizR)

#Rosenbrock function
objective <- function(x) {
  return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
}

#Rosenbrock gradient function
gradient <- function(x) {
  return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]), 200 * (x[2] - x[1] * x[1]) ) )
}

#estimated parameters
x<-c(-1.2, 0)


#minimizR function
opt<-minimizR(x,                   #initial parameters values
             objective,            #objective function
             gradient,             #gradient function
             control = list(       #control list
             tolerance = 1e-4,     #convergence criteria
             verbose = TRUE,       #print status
             iprint = 10,          #print interval
             hessian = TRUE))      #include hessian in the output 
```

## Output

```r
print(opt)
```

```
## $method
## [1] "l-bfgs"
## 
## $converged
## [1] TRUE
## 
## $`bounded problem`
## [1] FALSE
## 
## $message
## [1] "converged"
## 
## $iterations
## [1] 43
## 
## $`runtime (seconds)`
## [1] 0.030909
## 
## $`function value`
## [1] 2.166398e-11
## 
## $`norm gradient`
## [1] 7.353189e-05
## 
## $`max gradient component`
## [1] 6.744098e-05
## 
## $gradient
## [1] -6.744098e-05  2.930278e-05
## 
## $hessian
##           [,1]      [,2]
## [1,]  801.9929 -399.9982
## [2,] -399.9982  200.0000
## 
## $`det of hessian`
## [1] 399.9891
## 
## $`parameter values`
## [1] 0.9999956 0.9999913
```


