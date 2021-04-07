library(glmnet)
# lasso.included: TRUE  --> also returns lasso results
#                 FALSE --> does not return lasso results

# alpha.weights: 1 --> initial estimator : lasso --> onestep
#                0 --> initial estimator : ridge 
#               99 --> initial estimator : ols   

Adalasso.Nestedcv <- function(x, y, lasso.included=TRUE, alpha.weights=1, eps=0, nested.foldid=NULL, nfolds=10, family="gaussian"){
  # standard Lasso
  if (lasso.included==TRUE || alpha.weights==1) {
      Lasso =    Nestedcv.glmnet(nested.cv=FALSE, x=x, y=y, family=family, penalty.factor=rep(1,ncol(x)), eps=eps, alpha.weights=1,             nested.foldid=nested.foldid, nfolds=nfolds, std=TRUE , Index.nzero=1:ncol(x))
  }
  # Adaptive Lasso  
  if(alpha.weights==99){
      modstep1 = glm(y~x,family=family)
      Weights=1/(abs(modstep1$coefficients[-1])+eps)
      Index.nzero=which(modstep1$coefficients[-1]!=0)
  }else{
      if (alpha.weights==1) {
      modstep1 = Lasso  
      }else{
      modstep1 = Nestedcv.glmnet(nested.cv=FALSE, x=x, y=y, family=family, penalty.factor=rep(1,ncol(x)), eps=eps, alpha.weights=alpha.weights, nested.foldid=nested.foldid, nfolds=nfolds, std=TRUE , Index.nzero=1:ncol(x))
      }
      Weights=1/(abs( coef(modstep1, s= modstep1$lambda.min)[-1]) + eps)
      Index.nzero=which(coef(modstep1, s= modstep1$lambda.min)[-1]!=0)
   } 
      if (length(Index.nzero)!=0) {
      modstep2 = Nestedcv.glmnet(nested.cv=TRUE,  x=x, y=y, family=family, penalty.factor=Weights,        eps=eps, alpha.weights=alpha.weights, nested.foldid=nested.foldid, nfolds=nfolds, std=FALSE, Index.nzero=Index.nzero)
      }else{modstep2="Initial estimator returns null coefficients"}
  
    ifelse(lasso.included,return(list(Lasso=Lasso,Adaptive.Lasso=modstep2)),return(Adaptive.Lasso=modstep2))
}


# nested.cv : FALSE --> simple sheme
#             TRUE  --> nested sheme


Nestedcv.glmnet <- function(nested.cv=FALSE, x, y, family="gaussian", penalty.factor=rep(1,ncol(x)), eps=0, alpha.weights=1, nested.foldid=NULL, nfolds=10,std=TRUE,Index.nzero=1:ncol(x)){
  

  if(is.null(nested.foldid)){nested.foldid=rep_len(1:nfolds,length(y))} # for randomize: sample(rep_len(1:nfolds,length(y)))
  
  if (length(Index.nzero)==1) {
    Index.nzero=1:ncol(x)   
  }
  
  mod =     cv.glmnet(x=scale(x[,Index.nzero], center=FALSE, scale=penalty.factor[Index.nzero]), y=y, family=family, alpha = ifelse(nested.cv,1,alpha.weights), foldid=nested.foldid,standardize = std)
  
  if (nested.cv==TRUE) {
    
    nested.cverror=          matrix(NA,nrow = max(nested.foldid),ncol = length(mod$lambda))
  
    for(k in 1:max(nested.foldid)) {

      nested.fold <- which(nested.foldid == k)
      
      if(alpha.weights==99){
        mod.fold = glm(y[-nested.fold]~x[-nested.fold,],family=family)
        Weights=1/(abs(mod.fold$coefficients[-1])+eps)
        Index.nzero.fold=which(mod.fold$coefficients[-1]!=0)
      }else{
        mod.fold = cv.glmnet(x=x[-nested.fold,], y=y[-nested.fold], family=family, alpha = alpha.weights, penalty.factor=rep(1,ncol(x)),foldid = rep_len(1:nfolds,length(y[-nested.fold])))
        Weights= 1/(abs( coef(mod.fold, s= mod.fold$lambda.min)[-1]) + eps)
        Index.nzero.fold=which(coef(mod.fold, s= mod.fold$lambda.min)[-1]!=0)
      }
      
        
      if (length(Index.nzero.fold)==0) {
        pred = as.matrix(rep(mean(y[-nested.fold]),length(y[nested.fold])))
      }else{
        
      if (length(Index.nzero.fold)==1) {
        Index.nzero.fold=1:ncol(x[-nested.fold,])
      }
      
      modad.fold  = glmnet(x= scale(x[-nested.fold,Index.nzero.fold], center=FALSE, scale=Weights[Index.nzero.fold]),y=y[-nested.fold],family=family,alpha = 1,standardize = FALSE)
      pred = predict(modad.fold,  newx=scale(x[nested.fold,Index.nzero.fold], center=FALSE, scale=Weights[Index.nzero.fold]), type = "response",s = mod$lambda)
    
      }
      
      nested.cverror[k,]    <- apply((y[nested.fold] - pred  )^2, 2, mean)
      
    }
   
    cvall         = nested.cverror
    cvm           = apply(cvall, 2, mean)
    cvsd          = apply(cvall, 2, sd)/sqrt(max(nested.foldid))
    cvup          = cvm + cvsd
    cvlo          = cvm - cvsd
    lambda.min    = mod$lambda[which.min(cvm)]
    under.cvup1se = which(cvm<=cvup[which.min(cvm)]*1)
    min.nzero1se  = under.cvup1se[which.min(mod$nzero[under.cvup1se])]
    lambda.1se    = mod$lambda[min.nzero1se]
    
    add.to.mod=c("cvall","cvm","cvsd","cvup","cvlo","lambda.min","lambda.1se")
    
    for (add in add.to.mod) {eval(parse(text = paste0("mod$nested.",add,"=",add)))}
    
    mod$glmnet.fit$beta=apply(mod$glmnet.fit$beta, 2, function(b){
      r=rep(0,length(penalty.factor))
      r[Index.nzero]=b/penalty.factor[Index.nzero]
      return(r)})
    
  }

  return(mod)
  
}


