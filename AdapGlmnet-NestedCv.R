# adap: TRUE  --> lasso adap 
#       FALSE --> lasso std

# weights.alpha: 1 --> initial estimator : lasso --> onestep
#                0 --> initial estimator : ridge 
#               99 --> initial estimator : mco   



adap.glmnet <- function(adap=TRUE,weights.alpha=1,eps=0,nested.foldid=NULL,x,y,family="gaussian",alpha =1,penalty.factor=rep(1,ncol(x)), intercept = FALSE){
  
  if (adap==FALSE) {
    modstep1=nested.cvglmnet(nested.cv=adap,nested.foldid=nested.foldid,x=x,y=y,family=family,alpha =alpha,penalty.factor=penalty.factor, intercept = intercept)
    return(modstep1) 
  }else{
    
    if(weights.alpha==99){
      #add intercept option
      modstep1 = glm(y~x,family=family)
      Weights=1/(abs(modstep1$coefficients[-1])+eps)
      
      modstep2 = nested.cvglmnet(penalty.factor=Weights,nested.cv=adap,weights.alpha=weights.alpha,eps=eps,nested.foldid=nested.foldid,x=x,y=y,family=family,alpha=alpha, intercept = intercept)
      
    }else{
      
      modstep1 = nested.cvglmnet(penalty.factor=rep(1,ncol(x)),nested.cv=!adap,nested.foldid=nested.foldid,x=x,y=y,family=family,alpha=alpha, intercept = FALSE)
      Weights=1/(abs(modstep1$glmnet.fit$beta[,which(modstep1$glmnet.fit$lambda==modstep1$lambda.min)])+eps)
      
      modstep2 = nested.cvglmnet(penalty.factor=Weights,nested.cv=adap,weights.alpha=weights.alpha,eps=eps,nested.foldid=nested.foldid,x=x,y=y,family=family,alpha=alpha, intercept = intercept)
    }
    return(modstep2) 
  }
  
}


# nested.cv : FALSE --> cv.glmnet std
#             TRUE  --> cv.glmnet + nested.cv



nested.cvglmnet <- function(nested.cv=FALSE,weights.alpha=1,eps=0,nested.foldid=NULL,x,y,family,alpha =1,penalty.factor=rep(1,ncol(x)),intercept=TRUE,nfolds=10){
  
  
  #Standard
  #if(is.null(nested.foldid)){nested.foldid=rep_len(1:nfolds,length(y))[sample(length(y))]}
  
  #For Simulations
  if(is.null(nested.foldid)){nested.foldid=rep_len(1:nfolds,length(y))}
  
  
  mod = cv.glmnet(x=x,y=y,family=family,alpha = alpha,penalty.factor=penalty.factor,intercept=intercept,foldid=nested.foldid)
  
  if (nested.cv==TRUE) {
    
    
    nested.cverror=matrix(NA,nrow = max(nested.foldid),ncol = length(mod$lambda))
    
    
    
    for(k in 1:max(nested.foldid)) {
      
      nested.fold <- which(nested.foldid == k)
      
      if(weights.alpha==99){
        
        mod.fold = glm(y[-nested.fold]~x[-nested.fold,],family=family)
        Weights=1/(abs(mod.fold$coefficients[-1])+eps)
        
      }else{
        
        mod.fold = cv.glmnet(x=x[-nested.fold,],y=y[-nested.fold],family=family,alpha = weights.alpha,penalty.factor=rep(1,ncol(x)),intercept=intercept)
        Weights=1/(abs(mod.fold$glmnet.fit$beta[,which(mod.fold$glmnet.fit$lambda==mod.fold$lambda.min)])+eps)
      }
      
      modad.fold =glmnet(x=x[-nested.fold,],y=y[-nested.fold],family=family,lambda=mod$lambda,penalty.factor = Weights,alpha = alpha,intercept=intercept)
      
      #Missing: generalize the cv.error to all family (logistic, multinomial,...)
      for (l in 1:length(mod$lambda)) {
        nested.cverror[k,l]= (1/length(y[nested.fold]))*sum((y[nested.fold]-(x[nested.fold,]%*%modad.fold$beta[,l]+modad.fold$a0[l]))^2)
      }
      
    }
    
    cvm=apply(nested.cverror, 2, mean)
    cvsd=apply(nested.cverror, 2, sd)/sqrt(max(nested.foldid))
    cvup=cvm+cvsd
    cvlo=cvm-cvsd
    lambda.min=mod$lambda[which.min(cvm)]
    
    under.cvup1se=which(cvm<=cvup[which.min(cvm)]*1)
    min.nzero1se=under.cvup1se[which.min(mod$nzero[under.cvup1se])]
    
    lambda.1se=mod$lambda[min.nzero1se]
    
    addto.glmnet=c("cvm","cvsd","cvup","cvlo","lambda.min","lambda.1se")
    
    for (add in addto.glmnet) {eval(parse(text = paste0("mod$nested.",add,"=",add)))}
    
  }
  return(mod)
}


