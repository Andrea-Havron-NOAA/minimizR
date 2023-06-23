# minimizR
minimizR is a l-bfgs function minimizer that works with freely estimable or value bounded parameter sets.

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

## Value Bounded Example


Here is an example of minimizing the Rosenbrock function. Note, this example uses value bounded parameters. If "lb" and "ub" are omitted from the control list, all parameters will be freely estimated.  

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
opt<-minimizR(x,                    #initial parameters values
              objective,            #objective function
              gradient,             #gradient function
              control = list(       #control list
              tolerance = 1e-4,     #convergence criteria
              lb = c(-1.5, -1.5),   #lower bounds
              ub = c(1.5,1.5),      #upper bounds
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
## [1] TRUE
## 
## $message
## [1] "converged"
## 
## $iterations
## [1] 39
## 
## $`runtime (seconds)`
## [1] 0.01882
## 
## $`function value`
## [1] 1.075174e-14
## 
## $`norm gradient`
## [1] 3.058532e-06
## 
## $`max gradient component`
## [1] 2.767141e-06
## 
## $gradient
## [1] -2.767141e-06  1.302899e-06
## 
## $hessian
##           [,1] [,2]
## [1,]  801.9999 -400
## [2,] -400.0000  200
## 
## $`det of hessian`
## [1] 400.0003
## 
## $`parameter values`
## [1] 0.9999999 0.9999998
```


# Using minimizR with Template Model Builder

## About Template Model Builder

Template Model Builder is a framework that uses Automatic Differentiation to fit models to data. More information can be found at <https://cran.r-project.org/web/packages/TMB/index.html>. 

## R Code Using minimizR with TMB

Below is the TMB thetalog example using minimizR. The original code can be found at <https://github.com/kaskr/adcomp/tree/master/tmb_examples>.

```r
library(minimizR)
library(TMB)
compile("tmb_src/thetalog.cpp")
dyn.load(dynlib("tmb_src/thetalog"))

## Read data
Y <- scan("tmb_src/thetalog.dat", skip=3, quiet=TRUE)
data <- list(Y=Y)

## Parameter initial guess
parameters <- list(
  X = data$Y*0,
  logr0 = 0,
  logtheta = 0,
  logK = 6,
  logQ = 0,
  logR = 0
)

## Create the AD Function
obj <- MakeADFun(data, parameters, random="X", DLL="thetalog", silent = TRUE)
newtonOption(obj, smartsearch=FALSE)
```

## Initial Function Value And Gradient

```r
obj$fn()
```
```
## [1] 300.1637
## attr(,"logarithm")
## [1] TRUE
```
```r
obj$gr()
```
```
## [1]  42.40982  59.52237 -82.23018  40.40233  27.83362
```

## Run The Minimizer and View The Results

```r
#Fit the model
opt <- minimizR(obj$par, obj$fn, obj$gr, control = list(hessian = TRUE))
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
## [1] 16
## 
## $`runtime (seconds)`
## [1] 0.126925
## 
## $`function value`
## [1] 2.971193
## 
## $`norm gradient`
## [1] 1.545539e-05
## 
## $`max gradient component`
## [1] 1.493623e-05
## 
## $gradient
## [1]  3.263772e-06 -1.716953e-06  1.493623e-05  1.950193e-08  1.475689e-06
## 
## $hessian
##            [,1]        [,2]       [,3]      [,4]       [,5]
## [1,] 31.7912982   9.6840379  -1.825550 -6.408124  0.1369642
## [2,]  9.6840379   8.9235119 -14.165259 -4.414461 -0.9072788
## [3,] -1.8255537 -14.1652653 383.344716  4.909410  3.4215747
## [4,] -6.4081236  -4.4144607   4.909408 11.844641 10.5088531
## [5,]  0.1369642  -0.9072788   3.421574 10.508853 66.6376511
## 
## $`det of hessian`
## [1] 35511732
## 
## $`parameter values`
## [1] -2.6032942  0.7625685  6.7250075 -4.7496015 -3.1889238
```

## View The SD Report From TMB

```r
rep <- sdreport(obj)
print(rep)
```
```
## sdreport(.) result
##            Estimate Std. Error
## logr0    -2.6032942 0.22139111
## logtheta  0.7625685 0.44958978
## logK      6.7250075 0.05329932
## logQ     -4.7496015 0.35387920
## logR     -3.1889238 0.13401580
## Maximum gradient component: 1.493623e-05
```



