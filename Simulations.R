library(Matrix)
library(glmnet)
library(parallel)
path="C://Users/ballout/Desktop/CDESS/PAPIER-3/Codes/N/Codes/"

source(paste0(path,"AdapGlmnet-NestedCv.R"))

#settings

vsim=1:100
vp=c(100,500,1000)
vsignal=c(0.4,0.8,1.2,2.4)
vs0=c(0.05,0.1,0.25)
eps=1e-4
ntrain=c(1000,2000)
ntest=ntrain*10


code.meth=rbind(c(FALSE,rep(TRUE,3)),c(NA,1,0,99))
colnames(code.meth)=c("Lasso","Onestep","Adlasso.ridge","Adlasso.mco")
rownames(code.meth)=c("adap","weights.alpha")

for(n in 1:length(ntrain)){
  for(p in vp){
    for(s0 in vs0){
      for(signal in vsignal){
        funparall<- function(sim){
          res=list()
          
          ntot=ntrain[n]+ntest[n]
          x=matrix(rnorm(p*ntot),ncol = p)
          
          beta=rep(0,p)
          beta[sample(1:p,s0*p)]= (2*rbinom(n = s0*p,size = 1,prob = 0.5)-1)*signal

          y=x%*%beta+rnorm(ntot,0,1)
          
          xtrain=x[1:ntrain[n],]
          ytrain=y[1:ntrain[n]]
          
          xtest=x[(ntrain[n]+1):ntot,]
          ytest=y[(ntrain[n]+1):ntot]
          
          if(p>=(ntrain[n]/2)){vweights.alpha=c(1,0)}else{vweights.alpha=c(1,0,99)}
          
          
         for (adap in c(FALSE,TRUE)) {
            if (adap==TRUE){
            for (weights.alpha in vweights.alpha) {
            
             mod=adap.glmnet(adap=adap,weights.alpha=weights.alpha,eps=0,x=xtrain,y=ytrain,family="gaussian",alpha =1,intercept = FALSE)
             
             eval(parse(text=paste0("res$",names(which(code.meth[1,]==adap & code.meth[2,]==weights.alpha)),"=mod")))}
              }else{
              
              mod=adap.glmnet(adap=adap,x=xtrain,y=ytrain,family="gaussian",alpha =1,intercept = FALSE)
              
              res$Lasso=mod}
         }
          
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
