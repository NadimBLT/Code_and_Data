library(Matrix)
library(glmnet)
library(parallel)
library(here)

#path="C://Users/ballout/Desktop/CDESS/PAPIER-3/Codes/N/Codes/"

source("AdapGlmnet-NestedCv.R")

#settings

vsim    = 1:100
vp      = c(100,500,1000)
vsignal = c(0.25,0.5,1,2)
vs0     = c(10, 50)
eps     = 1e-4
ntrain  = c(1000,2000)
ntest   = 10000


code.meth=rbind(c(FALSE,rep(TRUE,3)),c(NA,1,0,99))
colnames(code.meth)=c("Lasso","Onestep","Adlasso.ridge","Adlasso.mco")
rownames(code.meth)=c("adap","alpha.weights")

for(n in 1:length(ntrain)){
  for(p in vp){
    for(s0 in vs0){
      for(signal in vsignal){
        
        # n= 1; p=100; s0=10; signal = 0.5
        funparall<- function(sim){
          res=list()
          
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
             
             eval(parse(text=paste0("res$",names(which(code.meth[1,]==adap & code.meth[2,]==alpha.weights)),"=mod")))}
              }else{
              
              mod=adap.glmnet(adap=adap,x=xtrain,y=ytrain,family="gaussian")#,alpha =1, intercept = TRUE)
              
              res$Lasso=mod}
         }
          
          # plot(res$Onestep$lambda[30:100], res$Onestep$cvm[30:100], log="x")
          # plot(res$Onestep$lambda[30:100], res$Onestep$nested.cvm[30:100], log="x")  
          # 
          # plot(res$Onestep$lambda, res$Onestep$cvm, log="x")
          # plot(res$Onestep$lambda, res$Onestep$nested.cvm, log="x")  
          # 
          # 
          # cbind(beta, res$Lasso$glmnet.fit$beta[, which(res$Lasso$lambda == res$Lasso$lambda.min)], res$Onestep$glmnet.fit$beta[, which(res$Onestep$lambda == res$Onestep$lambda.min)], 
          #       res$Onestep$glmnet.fit$beta[, which(res$Onestep$lambda == res$Onestep$nested.lambda.min)])
          # 
          
        res$train$x=xtrain
        res$train$y=ytrain
        
        res$test$x=xtest
        res$test$y=ytest
        
        res$beta.star=beta
        
        # we can do the evaluations here to avoid saving large files
        # ex: ntrain = 1000 + p = 20  ~ save = 1.5 mb 
        #                   + p = 100 ~ save = 8.0 mb
        #                   + p = 500 ~ save = 40.0 mb
         
        save(res,file = paste0(path,"res/ntrain-",ntrain[n],"-ntest-",ntest[n],"-nvar-",p,"-sparsity-",s0,"-signal-",signal,"sim-",sim,".RData"))
       
        
        }
        mclapply(vsim,FUN = funparall,mc.cores = 10)
        
      }
    }
  }
}
