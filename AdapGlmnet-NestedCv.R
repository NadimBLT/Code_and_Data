# adap: TRUE  --> lasso adap 
#       FALSE --> lasso std

# alpha.weights: 1 --> initial estimator : lasso --> onestep
#                0 --> initial estimator : ridge 
#               99 --> initial estimator : mco   



adap.glmnet <- function(x, y, adap=TRUE, alpha.weights=1, eps=0, nested.foldid=NULL, nfolds=10, family="gaussian", penalty.factor=rep(1,ncol(x))){
  
  if (adap==FALSE) {# standard Lasso
    modstep1=nested.cvglmnet(nested.cv=adap,nested.foldid=nested.foldid,x=x,y=y,family=family,penalty.factor=penalty.factor)
    return(modstep1) 
  }else{# Adaptive Lasso
    
  if(alpha.weights==99){ # OLS-AdaLasso
      modstep1 = glm(y~x,family=family)
      Weights=1/(abs(modstep1$coefficients[-1])+eps)
      
      # Selection of the parameter using (Nested and not std) CV
      modstep2 = nested.cvglmnet(penalty.factor=Weights,nested.cv=adap,alpha.weights=alpha.weights,eps=eps,nested.foldid=nested.foldid,x=x,y=y,family=family)
      
    }else{# weights derived from glmnet, with alpha set to alpha.weights (1: Lasso-> one-step Lasso; 0: ridge -> Ridge-AdaLasso)
      # in modstep1 below, you forgot to add the alpha.weights... so these weights were always computed with the Lasso, ...
      # that may explain the issues you observed for the Ridge-AdaLasso...
      
      #modstep1 = nested.cvglmnet(penalty.factor=rep(1,ncol(x)), nested.cv=!adap, nested.foldid=nested.foldid,x=x,y=y,family=family)
      #Weights=1/(abs(modstep1$glmnet.fit$beta[,which(modstep1$glmnet.fit$lambda==modstep1$lambda.min)])+eps)
      modstep1 = cv.glmnet(x=x,y=y,family=family, alpha = alpha.weights, penalty.factor=rep(1,ncol(x)), nfold=ifelse(is.null(nested.foldid), nfolds, max(nested.foldid)))
      Weights=1/(abs( coef(modstep1, s= modstep1$lambda.min)[-1]) + eps)
      
      
      
      modstep2 = nested.cvglmnet(penalty.factor=Weights,nested.cv=adap,alpha.weights=alpha.weights,eps=eps,nested.foldid=nested.foldid,x=x,y=y,family=family)
    }
    return(modstep2) 
  }
  
}


# nested.cv : FALSE --> cv.glmnet std
#             TRUE  --> cv.glmnet + nested.cv



nested.cvglmnet <- function(nested.cv=FALSE, alpha.weights=1, eps=0, nested.foldid=NULL, x, y, family="gaussian", penalty.factor=rep(1,ncol(x)), nfolds=10){
  
  # creation of the fold id for the cross-validation
  if(is.null(nested.foldid)){nested.foldid=sample(rep_len(1:nfolds,length(y)))}
  
  # first CV (the only one that is usually applied)
  mod = cv.glmnet(x=x, y=y, family=family, penalty.factor=penalty.factor, foldid=nested.foldid, alpha = 1)
  
  
  if (nested.cv==TRUE) {# second, "nested"-CV

    nested.cverror=matrix(NA,nrow = max(nested.foldid),ncol = length(mod$lambda))
    
    for(k in 1:max(nested.foldid)) {
      # we first re-compute the weights on each train sample
      # using either OLS fit (OLS-AdaLasso) or glmnet fit (one-step Lasso or Ridge-AdaLasso)

      nested.fold <- which(nested.foldid == k)
      
      if(alpha.weights==99){# OLS fit for the weights
        
        mod.fold = glm(y[-nested.fold]~x[-nested.fold,],family=family)
        Weights=1/(abs(mod.fold$coefficients[-1])+eps)
        
      }else{ # glmnet fit for the weights, with alpha set to alpha.weights
        
        mod.fold = cv.glmnet(x=x[-nested.fold,], y=y[-nested.fold], family=family, alpha = alpha.weights, penalty.factor=rep(1,ncol(x)))
        #Weights=1/(abs(mod.fold$glmnet.fit$beta[,which(mod.fold$glmnet.fit$lambda==mod.fold$lambda.min)])+eps)
        Weights= 1/(abs( coef(mod.fold, s= mod.fold$lambda.min)[-1]) + eps)
        
      }
      
      # compute the adaptive Lasso on the train sample, with penalty.factor set to Weights (computed on the train sample), 
      # and the lambda sequence computed from the whole sample
      modad.fold =glmnet(x=x[-nested.fold,],y=y[-nested.fold],family=family, lambda=mod$lambda, penalty.factor = Weights,alpha = 1)
      
      # And now compute the CV prediction error on the test sample
      #  => to be done: generalize the cv.error to all family (logistic, multinomial,...)
      for (l in 1:length(mod$lambda)) {
        nested.cverror[k,l]= (1/length(y[nested.fold]))*sum((y[nested.fold]-(x[nested.fold,]%*%modad.fold$beta[,l]+modad.fold$a0[l]))^2)
      }
      
    }
    
    cvm  = apply(nested.cverror, 2, mean)
    cvsd = apply(nested.cverror, 2, sd)/sqrt(max(nested.foldid))
    cvup = cvm + cvsd
    cvlo = cvm - cvsd
    lambda.min    = mod$lambda[which.min(cvm)]
    under.cvup1se = which(cvm<=cvup[which.min(cvm)]*1)
    min.nzero1se  = under.cvup1se[which.min(mod$nzero[under.cvup1se])]
    lambda.1se    = mod$lambda[min.nzero1se]
    
    addto.glmnet=c("cvm","cvsd","cvup","cvlo","lambda.min","lambda.1se")
    
    for (add in addto.glmnet) {eval(parse(text = paste0("mod$nested.",add,"=",add)))}
    
  }
  return(mod)
}


