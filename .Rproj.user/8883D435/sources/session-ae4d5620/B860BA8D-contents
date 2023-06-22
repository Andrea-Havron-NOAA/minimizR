# minimizR
minimizR is a l-bfgs function minimizer that works with freely estimable or value bounded parameter sets.

# Install from R

```
remotes::install_github("nmfs-fish-tools/minimizR")
library(minimizR)
```



## Freely Estimated Example

Here is an example of minimizing the Rosenbrock function with freely estimated parameters.

```{r min}
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

```{r out}
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
## [1] 0.038289
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
## $`parameter values`
## [1] 0.9999956 0.9999913
```

## Value Bounded Example


Here is an example of minimizing the Rosenbrock function. Note, this example uses value bounded parameters. If "lb" and "ub" are omitted from the control list, all parameters will be freely estimated.  

```{r minb}
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

```{r outb}
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
## [1] 0.013369
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
## $`parameter values`
## [1] 0.9999999 0.9999998
```
