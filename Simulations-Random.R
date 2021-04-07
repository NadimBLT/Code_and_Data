library(Matrix)
library(glmnet)
library(parallel)
library(here)

setwd("~")

source("Adalasso.Nestedcv.R")



Eval <- function(beta.star,beta.hat,intercept.hat,x,y,signal){
  

  Deviance   =  sum((y - (x%*%beta.hat + intercept.hat))^2)
  MSE        = mean((y - (x%*%beta.hat + intercept.hat))^2)
  
  EstErr1    =      sum((beta.hat-beta.star)^2)/sum((beta.star)^2)
  EstErr2    =      sum((beta.hat-beta.star)^2)/signal
  EstErr3    = sqrt(sum((beta.hat-beta.star)^2))
  
  beta.star[which(beta.star<0)]=-1
  beta.star[which(beta.star>0)]=1
  beta.hat[which(beta.hat<0)]=-1
  beta.hat[which(beta.hat>0)]=1
  
  Accuracysign = length(which(beta.hat==beta.star))/length(beta.star)
  Accuracypen  = length(which(beta.hat==beta.star))/length(beta.star) - length(which(beta.hat*beta.star==-1))/length(beta.star)
  
  beta.star[which(beta.star!=0)]=1
  beta.hat[which(beta.hat!=0)]=1
  
  Table=table(beta.star,beta.hat)
  if (ncol(Table)==1) {
    Table=cbind(Table,0) 
  }
  
  TP=Table[2,2]
  TN=Table[1,1]
  FP=Table[1,2]
  FN=Table[2,1]
  
  Accuracy  = (TP+TN)/(TP+TN+FP+FN)
  Precision = TP/(TP+FP)
  Recall    = TP/(TP+FN)
  F1score   = 2*(Recall*Precision)/(Recall+Precision)
  

  return(c(Accuracysign,Accuracypen,Accuracy,Precision,Recall,F1score,MSE,Deviance,EstErr1,EstErr2,EstErr3))
}


#settings
NCores  = 1 # number of cores (used for the parallel computation below)
vsim    = 1:100
vp      = c(100,500,1000)
vs0     = c(5,10,50)
vsignal = c(0.25,0.5,1,1.5)
veps    = 0
ntrain  = c(100,500,1000)
ntest   = 10000
rho=0.3

TableALL=NULL

for(eps in veps){
 for(n in 1:length(ntrain)){
   ntot=ntrain[n]+ntest
  for(p in vp){
    for(s0 in vs0){
      for(signal in vsignal){
    
        funparall<- function(sim){
          # to save all
          ## res=list()
          
          evaluation <- matrix(ncol = 12,nrow = 0)
          colnames(evaluation)=c("Accuracy-Sign","Accuracy-Sign-Pen","Accuracy","Precision","Recall","F1-Score","MSE","Deviance","Est-Err1","Est-Err2","Est-Err3","Method")
          
          #x=matrix(rnorm(p*ntot),ncol = p)
          
          x= lapply(1:ntot, FUN = function(i){arima.sim(p,model = list(ar=rho),sd=0.25)})
          #x= as.matrix(Reduce("rbind",x))
          x= matrix(unlist(x),ncol = p,byrow = TRUE)
              
          beta=rep(0,p)
          beta[sample(1:p,s0)]= (2*rbinom(n = s0,size = 1,prob = 0.5)-1)*signal

          y=x%*%beta+rnorm(ntot,0,1)
         
          xtrain=x[1:ntrain[n],]
          ytrain=y[1:ntrain[n]]
          
          xtest=x[(ntrain[n]+1):ntot,]
          ytest=y[(ntrain[n]+1):ntot]
          
     
          Adalasso.lasso.weights = Adalasso.Nestedcv(xtrain, ytrain, lasso.included=TRUE,  alpha.weights=1,  eps=0, nested.foldid=NULL, nfolds=10, family="gaussian")
          lasso=coef(Adalasso.lasso.weights$Lasso,s=Adalasso.lasso.weights$Lasso$lambda.min)
          if (length(which(lasso[-1]!=0))==0) {
          adalasso.lw.nai=adalasso.lw.nest=lasso
          }else{
          adalasso.lw.nai =coef(Adalasso.lasso.weights$Adaptive.Lasso,s=Adalasso.lasso.weights$Adaptive.Lasso$lambda.min)
          adalasso.lw.nest=coef(Adalasso.lasso.weights$Adaptive.Lasso,s=Adalasso.lasso.weights$Adaptive.Lasso$nested.lambda.min)
          }
          
          Adalasso.ridge.weights = Adalasso.Nestedcv(xtrain, ytrain, lasso.included=FALSE, alpha.weights=0,  eps=0, nested.foldid=NULL, nfolds=10, family="gaussian")
          adalasso.rw.nai =coef(Adalasso.ridge.weights,s=Adalasso.ridge.weights$lambda.min)
          adalasso.rw.nest=coef(Adalasso.ridge.weights,s=Adalasso.ridge.weights$nested.lambda.min)
          
          if (p<ntrain[n]) {
          Adalasso.ols.weights   = Adalasso.Nestedcv(xtrain, ytrain, lasso.included=FALSE, alpha.weights=99, eps=0, nested.foldid=NULL, nfolds=10, family="gaussian")
          adalasso.ow.nai =coef(Adalasso.ols.weights,s=Adalasso.ols.weights$lambda.min)
          adalasso.ow.nest=coef(Adalasso.ols.weights,s=Adalasso.ols.weights$nested.lambda.min)  
          }else{
          Adalasso.ols.weights=NULL
          }



             
          evaluation <- rbind(evaluation,
                                 c(Eval(beta.star = beta,beta.hat = lasso[-1],           intercept.hat = lasso[1]           ,y=ytest,x=xtest,signal),"lasso"),
                                 c(Eval(beta.star = beta,beta.hat = adalasso.lw.nai[-1], intercept.hat = adalasso.lw.nai[1], y=ytest,x=xtest,signal),"adalasso.lw.nai"),
                                 c(Eval(beta.star = beta,beta.hat = adalasso.lw.nest[-1],intercept.hat = adalasso.lw.nest[1],y=ytest,x=xtest,signal),"adalasso.lw.nest"),
                                 c(Eval(beta.star = beta,beta.hat = adalasso.rw.nai[-1], intercept.hat = adalasso.rw.nai[1], y=ytest,x=xtest,signal),"adalasso.rw.nai"),
                                 c(Eval(beta.star = beta,beta.hat = adalasso.rw.nest[-1],intercept.hat = adalasso.rw.nest[1],y=ytest,x=xtest,signal),"adalasso.rw.nest")
                                 )
          if (p<ntrain[n]) {
            evaluation <- rbind(evaluation,
                                c(Eval(beta.star = beta,beta.hat = adalasso.ow.nai[-1],  intercept.hat = adalasso.ow.nai[1],  y=ytest,x=xtest,signal),"adalasso.ow.nai"),
                                c(Eval(beta.star = beta,beta.hat = adalasso.ow.nest[-1], intercept.hat = adalasso.ow.nest[1], y=ytest,x=xtest,signal),"adalasso.ow.nest")
            )  
          }
          
        evaluation <- data.frame(evaluation,eps= eps, p=p,s0=s0,signal=signal,n=ntrain[n],sim=sim)
         
        onerunres=list(data=list(xtrain=xtrain,xtest=xtest,ytrain=ytrain,ytest=ytest,beta=beta),res=list(Lasso=Adalasso.lasso.weights$Lasso,Adalasso.lasso.weights=Adalasso.lasso.weights$Adaptive.Lasso,Adalasso.ridge.weights=Adalasso.ridge.weights,Adalasso.ols.weights=Adalasso.ols.weights))
        save(onerunres,file = paste0("Details/ntrain-",ntrain[n],"-ntest-",ntest,"-eps-", eps, "-nvar-",p,"-sparsity-",s0,"-signal-",signal,"-Sim-",sim,".RData"))
        
        return(evaluation)
        }
        
       Res=mclapply(vsim,FUN = funparall,mc.cores = NCores)
       
       Table <- NULL
       for(i in 1:length(vsim)){Table=rbind(Table,Res[[i]])}
        
       #back up: record a table for each combination of p*s0*n*signal
        save(Table,file = paste0("Evaluation/ntrain-",ntrain[n],"-ntest-",ntest,"-eps-", eps, "-nvar-",p,"-sparsity-",s0,"-signal-",signal,"AllSim.RData"))
       
       TableALL=rbind(TableALL,Table)
      }
    }
  }
save(TableALL,file = paste0("Evaluation/Table",ntrain[n],".Rdata"))  
}  
}
save(TableALL,file = "Evaluation/TableALL.Rdata")

