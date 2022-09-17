## -----------------------------------------------------------------------------
## Title: R-function ada.lasso.nested.cv()
## -----------------------------------------------------------------------------
## Authors: Nadim Ballout, Lola Etievant, Vivian Viallon
## -----------------------------------------------------------------------------
## Description:  This function returns the optimal adaptive lasso estimator,
##               with hyperparameter selected by the "nested" cross-validation
##               scheme, in the case of the ols-adaptive lasso, the 
##               ridge-adaptive lasso and the one-step lasso
## -----------------------------------------------------------------------------
## BiomJ article: "On the use of cross-validation for the calibration of
##                 the adaptive lasso"
## -----------------------------------------------------------------------------
## Algorithms in the article:
##  |-----|--------------|------------|------------------|---------------------|
##  | Alg | Initial Est. | CV  Scheme |    Func. Input   |    Func. Output     |
##  |-----|--------------|------------|------------------|---------------------|
##  |  3  |    lasso     |   nested   | alpha.weights=1  | ada.lasso.nested.cv |
##  |-----|--------------|------------|------------------|---------------------|
##  |  5  |     ols      |   nested   | alpha.weights=99 | ada.lasso.nested.cv |
##  |-----|--------------|------------|------------------|---------------------|
##  |  7  |    ridge     |   nested   | alpha.weights=0  | ada.lasso.nested.cv |
##  |--------------------------------------------------------------------------|
## -----------------------------------------------------------------------------
## Required Package: glmnet
## -----------------------------------------------------------------------------
## Usage: ada.lasso.nested.cv(
##                           x,
##                           y,
##                           alpha.weights = 1,
##                           eps=  0,
##                           nfolds = 10
##                           )
## -----------------------------------------------------------------------------
## Arguments:
##
##  x               input matrix, of dimension n x p, where n is the number of
##                  observations and p is the number of variables; each row is
##                  an observation vector.
##
##  y               quantitative response variable.
##
##  alpha.weights   the parameter that controls the penalty of the initial
##                  estimator of the adaptive lasso estimator with
##                  alpha.weights=99 for OLS (non-penalty), alpha.weights=1 for
##                  the lasso penalty, alpha.weights = 0 for the ridge penalty
##                  and alpha.weights is between 0 and 1 for the elasticnet 
##                  penalty; default is alpha.weights = 1
##
##  eps             the parameter that can guarantee that each component of the
##                  weight vector is finite, where the weights are defined as
##                  w_j = 1 / (|beta_j| + eps) for j in [p]; default is eps = 0
##
##  nfolds          number of folds in cross-validation; default is nfolds = 10
## -----------------------------------------------------------------------------
## Value: list of output values
##
## ada.lasso.nested.cv  the optimal adaptive lasso estimator with hyperparameter
##                      selected by the "nested" cross-validation scheme
##
## ada.lasso            list of glmnet function output values for the adaptive
##                      lasso results where the number of lambda values is 100
##                      (a0: intercept sequence of length 100, beta: p x 100 
##                      matrix of coefficients, etc..); for more details about
##                      these values please see Value section in glmnet {glmnet}
##                      in the R Documentation 
##
## nested.cv            list of the ingredients of the nested cross-validation
##                      fit (cvm: the mean cross-validated error - a vector of
##                      length 100, cvsd: estimate of standard error of cvm, 
##                      cvup: upper curve = cvm + cvsd, cvlo: lower curve =
##                      cvm - cvsd, lambda.min: value of lambda that gives
##                      minimum cvm, lambda.1se: largest value of lambda such
##                      that error is within 1 standard error of the minimum)
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

ada.lasso.nested.cv <- function(x, y, alpha.weights = 1, eps = 0, nfolds = 10) {
  
  
  ##  Step 1: computation of the initial estimator and the weights
  if (alpha.weights == 99) {
    mod.step1 <- glm(y ~ x, family = "gaussian")
    
    weights <- 1 / (abs(mod.step1$coefficients[-1]) + eps)
    index.nzero <- which(mod.step1$coefficients[-1]!=0)
  } else {
    mod.step1 <- cv.glmnet(x = x, y = y, family = "gaussian",
                           alpha = alpha.weights,
                           foldid = rep_len(1:nfolds, length(y)), 
                           standardize = TRUE)
    
    weights <- 1 / (abs( coef(mod.step1, s = mod.step1$lambda.min)[-1]) + eps)
    index.nzero <- which(coef(mod.step1, s = mod.step1$lambda.min)[-1] != 0)
  }
  
  ##  Step 2: computation of the final estimates
  if (length(index.nzero) == 0) {
    mod.step2 <- "Initial estimator returns null coefficients"
    nested.cv <- NULL
    ada.lasso.nested.cv <- NULL
  } else {
    ### computation of the adaptive estimator
    
    if (length(index.nzero) == 1) {
      index.nzero <- 1:ncol(x)   
    }
    
    mod.step2 <- glmnet(x = scale(x[ ,index.nzero], center = FALSE, 
                                  scale = weights[index.nzero]),
                        y = y, family = "gaussian", alpha = 1,
                        foldid = rep_len(1:nfolds, length(y)), 
                        standardize = FALSE)
    
    mod.step2$beta <- apply(mod.step2$beta, 2,
                            function(b){
                              r=rep(0,length(weights))
                              r[index.nzero] =
                                b / weights[index.nzero]
                              return(r)
                            }
    )
    
    ### computation of the nested cross validation error
    nested.foldid <- rep_len(1:nfolds, length(y))
    nested.cverror <- matrix(NA, nrow = nfolds, ncol = length(mod.step2$lambda))
    
    for (k in 1:nfolds) {
      nested.fold <- which(nested.foldid == k)
      #### update of the weights
      if (alpha.weights == 99) {
        mod.fold <- glm(y[-nested.fold] ~ x[-nested.fold, ],
                        family = "gaussian")
        weights.fold <- 1 / (abs(mod.fold$coefficients[-1]) + eps)
        index.nzero.fold <- which(mod.fold$coefficients[-1] != 0)
      } else {
        mod.fold <- cv.glmnet(x = x[-nested.fold, ], y = y[-nested.fold], 
                              family="gaussian", alpha = alpha.weights, 
                              foldid = rep_len(1:nfolds,
                                               length(y[-nested.fold])))
        weights.fold <- 1 / (abs(coef(mod.fold,
                                      s = mod.fold$lambda.min)[-1]) + eps)
        index.nzero.fold <- which(coef(mod.fold, 
                                       s = mod.fold$lambda.min)[-1] != 0)
      }
      
      #### computation of the estimation error
      if (length(index.nzero.fold) == 0) {
        pred <- as.matrix(rep(mean(y[-nested.fold]), length(y[nested.fold])))
      } else {
        if (length(index.nzero.fold) == 1) {
          index.nzero.fold <- 1:ncol(x[-nested.fold, ])
        }
        
        mod.ada.fold  <- glmnet(x = scale(x[-nested.fold, index.nzero.fold], 
                                          center=FALSE, 
                                          scale=weights.fold[index.nzero.fold])
                                , y = y[-nested.fold], family = "gaussian", 
                                alpha = 1, standardize = FALSE)
        
        pred <- predict(mod.ada.fold,  
                        newx = scale(x[nested.fold, index.nzero.fold], 
                                     center=FALSE,
                                     scale=weights.fold[index.nzero.fold]), 
                        type = "response", s = mod.step2$lambda)
      }
      nested.cverror[k, ] <- apply((y[nested.fold] - pred  ) ^ 2, 2, mean)
    }
    
    ### calculation of the ingredients of the nested cross validation
    cvall         <- nested.cverror
    cvm           <- apply(cvall, 2, mean)
    cvsd          <- apply(cvall, 2, sd) / sqrt(max(nested.foldid))
    cvup          <- cvm + cvsd
    cvlo          <- cvm - cvsd
    lambda.min    <- mod.step2$lambda[which.min(cvm)]
    under.cvup1se <- which(cvm <= cvup[which.min(cvm)] * 1)
    min.nzero1se  <- under.cvup1se[which.min(mod.step2$nzero[under.cvup1se])]
    lambda.1se    <- mod.step2$lambda[min.nzero1se]
    
    nested.cv <- list()
    add.to.nested.cv <- c("cvm", "cvsd", "cvup", "cvlo", "lambda.min", 
                          "lambda.1se")
    for (add in add.to.nested.cv) {
      eval(parse(text = paste0("nested.cv$", add, " <- ", add)))}
    
    
    ada.lasso.nested.cv <- coef(mod.step2, s = nested.cv$lambda.min)
    
  }
  
  return(list(ada.lasso.nested.cv = ada.lasso.nested.cv,
              ada.lasso = mod.step2,
              nested.cv = nested.cv)
  )
}



