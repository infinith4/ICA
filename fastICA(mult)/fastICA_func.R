fastICAmult<-function(X,m,W=signaldata("unif",m,m),maxcnt=10000,epsilon=0.0001){
    if(dim(X)[2]>dim(X)[1]){ #行:siganl,列:標本
        X<-t(X);
    }
    
    V<-svd(var(X))$u%*%diag(1/sqrt(svd(var(X))$d))%*%t(svd(var(X))$u);
    Z<-V%*%t(X); #whiting
    #Z=V%*%X #ok
    #W <- matrix(scan("W_init01.txt"), nrow=m) #ok
    W<-t(W)
    #ノルム1にする
    for(i in 1:m){
        W[,i]<-W[,i]/vecnorm(W[,i])
    }
    i=1;
    while(i<=m){
        cat("Independent Component:",i,"\n");
        cnt=0;
        while(cnt<maxcnt){
            #cat("cnt:",cnt,"\n")
            #5.
            wbefore<-W[,i]
            W[,i]<-apply(t(t(Z)*tanh(t(W[,i])%*%Z)[1,]),1,mean)-mean(1-(tanh(t(W[,i])%*%Z)[1,])^2)*W[,i];
            #6.グラムシュミットの直交化
            if(2<=i){
                Sum=0;
                for(j in 1:(i-1)){
                    Sum<-Sum+(t(W[,i])%*%W[,j])*W[,j];
                }
                W[,i]<-W[,i]-Sum; #2<=i
            }
            #7.収束したか?
            W[,i]<-W[,i]/vecnorm(W[,i]);
            #cat("W[,i]:",W[,i],"\n\n")            
            if(1-epsilon<=abs(t(wbefore)%*%W[,i]) && abs(t(wbefore)%*%W[,i])<=1+epsilon){
                cat("cnt:",cnt,"\n");
                cat("W[,i]:",W[,i],"\n\n");
                i<-i+1;
                cnt=maxcnt;
            }
            cnt=cnt+1;
        }
    }
    
    return(list(V=V,W=W,Z=Z,restoredsignal=as.matrix(t(W)%*%Z[,1:n])))
}
