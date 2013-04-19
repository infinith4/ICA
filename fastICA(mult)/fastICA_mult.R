#複数の独立成分を推定するfastICA

library(splus2R);

original_unifdata <-function(p,n){
    S<-matrix(,p,1)
    i=0
    while(i<n){
        u<-runif(p,-sqrt(3),sqrt(3))
        S<-cbind(S,u)
        i<-i+1
    }
    S<-S[,-1]
    return(S)
}

n=10
p=5
S<-original_unifdata(p,n) #p,n
S
A=matrix(runif(p^2,-sqrt(3),sqrt(3)),p) #p,p
A
X=A%*%S #p,n

X
#sigma2 = mean((X[i,]-mean(X[i,]))^2)
#sqrt(sigma2)
#sd(X[1,])
#z=(X[i,]-mean(X[i,]))/sd(X[i,])
#mean(z)
#sd(z)

#z
#X=scale(t(X))[,]

#2.whiting
S<-original_unifdata(p,n) #p,n
X=A%*%S #p,n

X<-X-apply(X,2,mean)
X<-t(X)
X
var(X) #p,p
E<-eigen(var(X))$vectors
E
eigen(var(X))$values

sqrt(eigen(var(X))$values)
Dm12<-diag(1/sqrt(eigen(var(X))$values))
Dm12
V<-E%*%Dm12%*%t(E)
V%*%X[1,]
Z=V%*%t(X)

###################


m=5 #given
W<-diag(rep(1,m)) #初期値:単位行列 #given
Z #given

i=1
cnt=0
while(cnt<1000){
    #5.
    cat("cnt:",cnt,"\n")
    while(i<=m){
        W[,i]<-apply(t(t(Z)*tanh(t(W[,i])%*%Z)[1,]),1,mean)-mean(1-(tanh(t(W[,1])%*%Z)[1,])^2)*W[,i]
        i<-i+1
    }
    #W<-t(W)
    #6.
    E<-eigen(W%*%t(W))$vectors
    E=t(E)
    Dm12<-diag(sqrt(eigen(W%*%t(W))$values))
    epsilon=0.001
    if(1-epsilon<=vecnorm((E%*%Dm12%*%t(E)%*%W)[,1])&& vecnorm((E%*%Dm12%*%t(E)%*%W)[,1])<=1+epsilon){
        print("ok\n")
    }
    W<-E%*%Dm12%*%t(E)%*%W
    cat(W%*%t(W),"\n")
    cnt=cnt+1
}


apply(Z,1,sum) #Zの各行についての和

t(W[,1])%*%Z
Z

t(t(Z)*tanh(t(W[,1])%*%Z)[1,]) #ok:行が次元,列が要素数
apply(t(t(Z)*tanh(t(W[,1])%*%Z)[1,]),1,mean) #各行について平均

mean(1-(tanh(t(W[,1])%*%Z)[1,])^2)

apply(1-(tanh(t(W[,1])%*%Z)[1,])^2),1,mean)

Z
A

b=c(1:n^2)
b
B=matrix(b,p,n) #p,p
B
t(t(B)*c(1:10))



while(i<m){
    
    
}

fastICAs <- function(X,W_init,m){ #m:独立成分の数
    #1.標準化
    X=t(scale(t(X))[,])
    #2.白色化
    eigen(X)
    
}

fastICAs(n)