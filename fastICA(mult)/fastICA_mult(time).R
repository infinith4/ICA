#複数の独立成分を推定するfastICA

library(splus2R);
library(VGAM)
setwd("/home/th4/Dropbox/seminer/fastICA/fastICA(mult)")

source("signaldata.R")


n=500
p=3
S=signaldata("laplace",p,n)
plot(S[1,],xlim=c(0,n),ylim=c(min(S[1,]),max(S[1,])),type="l")
hist(S[1,])
S=signaldata("unif",p,n)
plot(S[1,],xlim=c(0,n),ylim=c(min(S[1,]),max(S[1,])),type="l")
hist(S[1,])

dim(signaldata("laplace",p,n))

#S<-signaldata("laplace",p,n) #p,n
#dim(S)
#write(S, file="./laplacedata01.txt")

#S<-signaldata("laplace",p,n) #p,n
#dim(S)
#write(S, file="./uniformdata01.txt")

# scan 関数で再読み込み(結果はベクトル)
#S <- matrix(scan("uniformdata01.txt"), ncol=n) #ok
S <- matrix(scan("laplacedata01.txt"), ncol=n) #ok
S[,1:2] #ok
dim(S)


#A<-matrix(runif(p^2,-sqrt(3),sqrt(3)),p)
#write(A, file="./A01.txt")
# scan 関数で再読み込み(結果はベクトル)
#A <- matrix(scan("A01.txt"), ncol=p)
#A<-t(A)

#A<-matrix(runif(n^2,-sqrt(3),sqrt(3)),n)
#dim(A)
#write(A, file="./A01.txt")
# scan 関数で再読み込み(結果はベクトル)
#A <- matrix(scan("A01.txt"), ncol=p)
A <- matrix(scan("A01.txt"), ncol=n)
A<-t(A)
dim(A)
dim(S)
S<-t(S)
#X=A%*%S #p,n
X=A%*%S #p,n
dim(X)
X[,1:2] #ok
X<-X-apply(t(X),2,mean) #
X<-t(X)

eigen(var(X))

#E<-eigen(var(X))$vectors
E<-eigen(var(X))$vectors[,1:3]
dim(E)
sqrt(eigen(var(X))$value[1:3])

#Dm12<-diag(1/sqrt(eigen(var(X))$values))
Dm12<-diag(1/sqrt(eigen(var(X))$values[1:p]))
dim(Dm12)
V<-E%*%Dm12%*%t(E) #ok
dim(V)

#2.whiting

#S<-original_unifdata(p,n) #p,n
#write(S, file="./whitinguniformdata01.txt")
#S<-signaldata("laplace",p,n) #p,n
#dim(S)
#write(S, file="./whitinglaplacedata01.txt")
# scan 関数で再読み込み(結果はベクトル)
#S <- matrix(scan("whitinguniformdata01.txt"), ncol=n)
S <- matrix(scan("whitinglaplacedata01.txt"), ncol=n)
S<-t(S)
X=A%*%S #p,n
X[,1:2]
X<-X-apply(t(X),2,mean) #ok
X[,1:2]
Z=V%*%X #ok
Z[,1:2]
dim(Z)

###################

m=3 #独立成分の数
W<-diag(rep(1,n)) #初期値:単位行列
#W<-signaldata("uniform",p,p) #p,n
#W<-original_unifdata(m,m)
#write(W, file="./W_init01.txt")
# scan 関数で再読み込み(結果はベクトル)
#W <- matrix(scan("W_init01.txt"), nrow=m) #ok
dim(W)
W<-t(W)

for(i in 1:m){
    W[,i]<-W[,i]/vecnorm(W[,i])
}
W
epsilon=0.0001
i=1
maxcnt=10000
while(i<=m){
    cat("i:",i,"\n")
    cnt=0
    while(cnt<maxcnt){
        #cat("cnt:",cnt,"\n")
        #5.
        wbefore<-W[,i]
        W[,i]<-apply(t(t(Z)*tanh(t(W[,i])%*%Z)[1,]),1,mean)-mean(1-(tanh(t(W[,i])%*%Z)[1,])^2)*W[,i]
        #6.グラムシュミットの直交化
        if(2<=i){
            Sum=0
            for(j in 1:(i-1)){
                Sum<-Sum+(t(W[,i])%*%W[,j])*W[,j]
            }
            W[,i]<-W[,i]-Sum #2<=i
        }
        #7.
        W[,i]<-W[,i]/vecnorm(W[,i])
        #cat("W[,i]:",W[,i],"\n\n")
        
        if(1-epsilon<=abs(t(wbefore)%*%W[,i]) && abs(t(wbefore)%*%W[,i])<=1+epsilon){
            cat("cnt:",cnt,"\n")
            cat("W[,i]:",W[,i],"\n\n")
            i<-i+1
            cnt=maxcnt
        }
        cnt=cnt+1
    }
}

Iresult <- W%*%t(W)
result <- t(W)%*%V%*%A
dim(result)
plot(result[1,],xlim=c(0,n),ylim=c(min(result[1,]),max(result[1,])),type="l")
plot(result[2,],xlim=c(0,n),ylim=c(min(result[2,]),max(result[2,])),type="l")
plot(result[3,],xlim=c(0,n),ylim=c(min(result[3,]),max(result[3,])),type="l")
hist(result[1,])
hist(result[2,])
hist(result[3,])

plot(S[,1],xlim=c(0,n),ylim=c(min(S[,1]),max(S[,1])),type="l")
hist(S[,1])
plot(S[,2],xlim=c(0,n),ylim=c(min(S[,2]),max(S[,2])),type="l")
hist(S[,2])
plot(S[,3],xlim=c(0,n),ylim=c(min(S[,3]),max(S[,3])),type="l")
hist(S[,3])
#################

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