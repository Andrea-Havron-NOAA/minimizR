library(minimizR)


Cos<-function(x){
    return(cos(x[1]))
}

Cos_dx<-function(x){
    return(-sin(x[1]))
}
x<-c(1.5)
opt<-minimize(x,Cos, Cos_dx)
print(opt)
