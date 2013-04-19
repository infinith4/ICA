#複数の独立成分を推定するfastICA

library(splus2R);
library(VGAM);
setwd("/home/th4/Dropbox/seminer/fastICA/fastICA(mult)");

source("signaldata.R");
source("fastICA_func.R");


n=10000
p=4
#dim(signaldata("unif",p,n))
S1<-signaldata("laplace",2,n)
S2<-signaldata("unif",2,n)
S<-rbind(S1,S2)
#dim(S)

plot(S[1,],type="l")
hist(S[1,],breaks="Scott",freq = FALSE)
#plot(density(X[1,]))
lines(density(S[1,]), col = "orange", lwd = 2)
rug(S[1,])

#  plot(S[2,],type="l")
#  hist(S[2,])
#  plot(S[3,],type="l")
#  hist(S[3,])

write(S, file="./laplacedata01.txt")
# scan 関数で再読み込み(結果はベクトル)
#S <- matrix(scan("uniformdata01.txt"), ncol=n) #ok
#S[,1:2] #ok
#dim(S)

A<-matrix(runif(p^2,-sqrt(3),sqrt(3)),p)
#write(A, file="./A01.txt")
# scan 関数で再読み込み(結果はベクトル)
#A <- matrix(scan("A01.txt"), ncol=p)
A<-t(A)
dim(A)

X=A%*%S #p,n
dim(X)
X[,1:2] #ok
X<-X-apply(t(X),2,mean) #
plot(X[1,],type="l")
hist(X[1,],breaks="Scott",freq = FALSE)
#plot(density(X[1,]))
lines(density(X[1,]), col = "orange", lwd = 2)
rug(X[1,])

plot(X[2,],type="l")
hist(X[2,])
plot(X[3,],type="l")
hist(X[3,])
X<-t(X)



###################

m=4 #独立成分の数
#W<-diag(rep(1,m)) #初期値:単位行列
#W<-original_unifdata(m,m)
#write(W, file="./W_init01.txt")
# scan 関数で再読み込み(結果はベクトル)


result<-fastICAmult(X,m)
V<-result$V
dim(V)

W<-result$W
W

Z<-result$Z

W%*%t(W)
t(W)%*%V%*%A

restoredsignal<-result$restoredsignal

#split.screen(c(1,3),screen=2)

par(mfrow=c(2,2))
#plot(restoredsignal[1,],type="l") #
for(i in 1:4){
    hist(restoredsignal[i,],breaks="Scott",freq = FALSE)
    #plot(density(X[1,]))
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