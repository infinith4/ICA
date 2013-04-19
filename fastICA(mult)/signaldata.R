library(VGAM)
signaldata <-function(func,p,n){
    S<-matrix(,p,1)
    i=0
    while(i<n){
        if(func == "unif"){
            u<-runif(p,-sqrt(3),sqrt(3))
        }else if(func == "laplace"){
            u<-rlaplace(p, location=0, scale=1)
        }
        S<-cbind(S,u)
        i<-i+1
    }
    S<-S[,-1]
    return(S)
}

