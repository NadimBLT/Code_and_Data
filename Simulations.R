library(Matrix)
library(glmnet)
library(parallel)
library(here)

#setwd("~")

source("AdapGlmnet-NestedCv.R")


Eval <- function(beta.star,beta.hat,intercept.hat,x,y){
  
  prediction.error= mean((y - (x%*%beta.hat + intercept.hat))^2)
  
  beta.star[which(beta.star<0)]=-1
  beta.star[which(beta.star>0)]=1
  beta.hat[which(beta.hat<0)]=-1
  beta.hat[which(beta.hat>0)]=1
  
  accuracy=length(which(beta.hat==beta.star))/length(beta.star)
  return(c(accuracy,prediction.error))
}

#settings

vsim    = 1:50
vp      = c(100,500,1000)
vsignal = c(0.25,0.5,1,1.5)
vs0     = c(10, 50)
eps     = 1e-4
ntrain  = c(1000)
ntest   = ntrain*10


code.meth=rbind(c(FALSE,rep(TRUE,3)),c(NA,1,0,99))
colnames(code.meth)=c("Lasso","One.step.Lasso","Ridge.AdaLasso","OLS.AdaLasso")
rownames(code.meth)=c("adap","alpha.weights")

TableALL=NULL

for(n in 1:length(ntrain)){
  for(p in vp){
    for(s0 in vs0){
      for(signal in vsignal){
    
        funparall<- function(sim){
          # to save all
          ## res=list()
          
          evaluation <- matrix(ncol = 4,nrow = 0)
          colnames(evaluation)=c("Accuracy","PredError","Method","Type")
          
          ntot=ntrain[n]+ntest[n]
          x=matrix(rnorm(p*ntot),ncol = p)
          
          beta=rep(0,p)
          beta[sample(1:p,s0)]= (2*rbinom(n = s0,size = 1,prob = 0.5)-1)*signal

          y=x%*%beta+rnorm(ntot,0,1)
          
          xtrain=x[1:ntrain[n],]
          ytrain=y[1:ntrain[n]]
          
          xtest=x[(ntrain[n]+1):ntot,]
          ytest=y[(ntrain[n]+1):ntot]
          
          if(p>=(ntrain[n]/2)){valpha.weights=c(1,0)}else{valpha.weights=c(1,0,99)}
          
          
         for (adap in c(FALSE,TRUE)) {
            if (adap==TRUE){
            for (alpha.weights in valpha.weights) {
            
             mod=adap.glmnet(adap=adap,alpha.weights=alpha.weights,eps=eps,x=xtrain,y=ytrain,family="gaussian")#,alpha =1, intercept = TRUE)
             
             # to save all :
             ## eval(parse(text=paste0("res$",names(which(code.meth[1,]==adap & code.meth[2,]==alpha.weights)),"=mod")))
            
             beta.hat1=coef(mod, s= mod$lambda.min)
             beta.hat2=coef(mod, s= mod$nested.lambda.min)
             
             evaluation <- rbind(evaluation,
                                 c(Eval(beta.star = beta,beta.hat = beta.hat1[-1],intercept.hat = beta.hat1[1],y=ytest,x=xtest),names(which(code.meth[1,]==adap & code.meth[2,]==alpha.weights)),"CV"),
                                 c(Eval(beta.star = beta,beta.hat = beta.hat2[-1],intercept.hat = beta.hat2[1],y=ytest,x=xtest),names(which(code.meth[1,]==adap & code.meth[2,]==alpha.weights)),"nestedCV"))
              
              }
              
              }else{
              
              mod=adap.glmnet(adap=adap,x=xtrain,y=ytrain,family="gaussian")#,alpha =1, intercept = TRUE)
              
              # to save all
              ## res$Lasso=mod
              
              beta.hat=coef(mod, s= mod$lambda.min)
              evaluation <- rbind(evaluation,c(Eval(beta.star = beta,beta.hat = beta.hat[-1],intercept.hat = beta.hat[1],y=ytest,x=xtest),"Lasso","CV"))
              
              }
         }
          
        evaluation <- data.frame(evaluation,p=p,s0=s0,signal=signal,n=ntrain[n],sim=sim)
         
        # to save all : 
        
        ##res$train$x=xtrain
        ##res$train$y=ytrain
        
        ##res$test$x=xtest
        ##res$test$y=ytest
        
        ##res$beta.star=beta
        
        ## res$evaluation <- evaluation
       
        ## save(res,file = paste0(path,"res/ntrain-",ntrain[n],"-ntest-",ntest[n],"-nvar-",p,"-sparsity-",s0,"-signal-",signal,"sim-",sim,".RData"))
        
        return(evaluation)
        
        }
        
       Res=mclapply(vsim,FUN = funparall,mc.cores = 50)
       
       Table <- NULL
       for(i in 1:length(vsim)){Table=rbind(Table,Res[[i]])}
        
       #back up: record a table for each combination of p*s0*n*signal
        save(Table,file = paste0("RES/ntrain-",ntrain[n],"-ntest-",ntest[n],"-nvar-",p,"-sparsity-",s0,"-signal-",signal,"AllSim.RData"))
       
       TableALL=rbind(TableALL,Table)
      }
    }
  }
}

save(TableALL,file = "RES/TableALL.Rdata")

