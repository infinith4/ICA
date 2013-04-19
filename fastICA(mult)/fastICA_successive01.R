#複数の独立成分を推定するfastICA

library(splus2R);
library(VGAM);
setwd("/home/th4/Dropbox/seminer/fastICA/fastICA(mult)");

source("signaldata.R");
source("fastICA_func.R");


n=10000
p=4
S1<-signaldata("laplace",2,n)
S2<-signaldata("unif",2,n)
S<-rbind(S1,S2)


par(mfrow=c(2,4))
for(i in 1:4){
    plot(S[i,],type="l")
    hist(S[i,],breaks="Scott",freq = FALSE)
    lines(density(S[i,]), col = "orange", lwd = 2)
    rug(S[i,])
}
A<-matrix(runif(p^2,-sqrt(3),sqrt(3)),p)
A<-t(A)
dim(A)

X=A%*%S #p,n
dim(X)
X<-X-apply(t(X),2,mean) #

par(mfrow=c(2,4))
for(i in 1:4){
    plot(X[i,],type="l")
    hist(X[i,],breaks="Scott",freq = FALSE)
    lines(density(X[i,]), col = "orange", lwd = 2)
    rug(X[i,])
}
X<-t(X)



###################

m=20 #独立成分の数

result<-fastICAmult(X,m)
V<-result$V
dim(V)

W<-result$W
W

Z<-result$Z
par(mfrow=c(2,4))
for(i in 1:4){
    plot(Z[i,],type="l")
    hist(Z[i,],breaks="Scott",freq = FALSE)
    lines(density(Z[i,]), col = "orange", lwd = 2)
    rug(Z[i,])
}

W%*%t(W)
t(W)%*%V%*%A #順序行列になれば成功

restoredsignal<-result$restoredsignal

par(mfrow=c(2,2))
for(i in 1:4){
    hist(restoredsignal[i,],breaks="Scott",freq = FALSE)
    lines(density(restoredsignal[i,]), col = "orange", lwd = 2)
    rug(restoredsignal[i,])
}


###################################################
S=0
i=2
for(j in 1:(m-1)){
    S<-S+(t(W[,i])%*%W[,j])*W[,j]
}
S

i=2;j=3
(t(W[,i])%*%W[,j])*W[,j]
W
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
library(fastICA)
fastICA(X,5)

fastICAs <- function(X,W_init,m){ #m:独立成分の数
    #1.標準化
    X=t(scale(t(X))[,])
    #2.白色化
    eigen(X)
    
}

fastICAs(n)