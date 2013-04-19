#複数の独立成分を推定するfastICA

library(splus2R);
library(VGAM)
setwd("/home/th4/Dropbox/seminer/fastICA/fastICA(mult)")

source("signaldata.R")


n=10000
p=5
dim(signaldata("unif",p,n))
S<-signaldata("unif"p,n) #p,n
#plot(S[1,],type="l")
#hist(S[1,])
#plot(S[2,],type="l")
#hist(S[2,])
#plot(S[3,],type="l")
#hist(S[3,])

#S<-signaldata(p,n) #p,n
#dim(S)
#write(S, file="./uniformdata01.txt")
# scan 関数で再読み込み(結果はベクトル)
S <- matrix(scan("uniformdata01.txt"), ncol=n) #ok
#S[,1:2] #ok
#dim(S)

#A<-matrix(runif(p^2,-sqrt(3),sqrt(3)),p)
#write(A, file="./A01.txt")
# scan 関数で再読み込み(結果はベクトル)
A <- matrix(scan("A01.txt"), ncol=p)
A<-t(A)
#dim(A)

X=A%*%S #p,n
#dim(X)
#X[,1:2] #ok
X<-X-apply(t(X),2,mean) #
X<-t(X)

#var(X)
E<-eigen(var(X))$vectors
Dm12<-diag(1/sqrt(eigen(var(X))$values))
V<-E%*%Dm12%*%t(E) #ok


#2.whiting

#S<-original_unifdata(p,n) #p,n
#write(S, file="./whitinguniformdata01.txt")
# scan 関数で再読み込み(結果はベクトル)
#S <- matrix(scan("whitinguniformdata01.txt"), ncol=n)
#X=A%*%S #p,n
#X[,1:2]
#X<-X-apply(t(X),2,mean) #ok
#X[,1:2]
Z=V%*%X #ok
Z[,1:2]
dim(Z)

###################

m=5 #独立成分の数
#W<-diag(rep(1,m)) #初期値:単位行列
#W<-original_unifdata(m,m)
#write(W, file="./W_init01.txt")
# scan 関数で再読み込み(結果はベクトル)
W <- matrix(scan("W_init01.txt"), nrow=m) #ok
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

W%*%t(W)
t(W)%*%V%*%A


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