## -----------------------------------------------------------------------------
## Title: R-function ada.lasso.naive.cv()
## -----------------------------------------------------------------------------
## Authors: Nadim Ballout, Lola Etievant, Vivian Viallon
## -----------------------------------------------------------------------------
## Description:  This function returns the optimal adaptive lasso estimator,
##               with hyperparameter selected by the "naive" cross-validation
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
##  |  2  |    lasso     |   naive    | alpha.weights=1  | ada.lasso.naive.cv  |
##  |-----|--------------|------------|------------------|---------------------|
##  |  4  |     ols      |   naive    | alpha.weights=99 | ada.lasso.naive.cv  |
##  |-----|--------------|------------|------------------|---------------------|
##  |  6  |    ridge     |   naive    | alpha.weights=0  | ada.lasso.naive.cv  |
##  |--------------------------------------------------------------------------|
## -----------------------------------------------------------------------------
## Required Package: glmnet
## -----------------------------------------------------------------------------
## Usage: ada.lasso.naive.cv(
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
## ada.lasso.naive.cv   the optimal adaptive lasso estimator with hyperparameter
##                      selected by the "naive" cross-validation scheme
##
## ada.lasso            list of glmnet function output values for the adaptive
##                      lasso results where the number of lambda values is 100
##                      (a0: intercept sequence of length 100, beta: p x 100 
##                      matrix of coefficients, etc..); for more details about
##                      these values please see Value section in glmnet {glmnet}
##                      in the R Documentation 
##
## naive.cv             list of the ingredients of the naive cross-validation
##                      fit (cvm: the mean cross-validated error - a vector of
##                      length 100, cvsd: estimate of standard error of cvm, 
##                      cvup: upper curve = cvm + cvsd, cvlo: lower curve =
##                      cvm - cvsd, lambda.min: value of lambda that gives
##                      minimum cvm, lambda.1se: largest value of lambda such
##                      that error is within 1 standard error of the minimum)
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

ada.lasso.naive.cv <- function(x, y, alpha.weights = 1, eps = 0, nfolds = 10) {
  
  
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
    naive.cv <- NULL
    ada.lasso.naive.cv <- NULL
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
                            function(b) {
                              r <- rep(0, length(weights))
                              r[index.nzero] <- b / weights[index.nzero]
                              return(r)
                            })
    
    ### computation of the naive cross validation error
    naive.foldid <- rep_len(1:nfolds, length(y))
    naive.cverror <- matrix(NA, nrow = nfolds, ncol = length(mod.step2$lambda))
    
    for (k in 1:nfolds) {
      naive.fold <- which(naive.foldid == k)
      
      #### computation of the estimation error
      
      if (length(index.nzero) == 1) {
        index.nzero <- 1:ncol(x[-naive.fold, ])
      }
      
      mod.ada.fold  <- glmnet(x = scale(x[-naive.fold, index.nzero], 
                                        center=FALSE, 
                                        scale=weights[index.nzero]),
                              y = y[-naive.fold], family = "gaussian", 
                              alpha = 1, standardize = FALSE)
      
      pred <- predict(mod.ada.fold,  
                      newx = scale(x[naive.fold, index.nzero], 
                                   center=FALSE,
                                   scale=weights[index.nzero]), 
                      type = "response", s = mod.step2$lambda)
      
      naive.cverror[k, ] <- apply((y[naive.fold] - pred  ) ^ 2, 2, mean)
    }
    
    ### calculation of the ingredients of the naive cross validation
    cvall         <- naive.cverror
    cvm           <- apply(cvall, 2, mean)
    cvsd          <- apply(cvall, 2, sd) / sqrt(max(naive.foldid))
    cvup          <- cvm + cvsd
    cvlo          <- cvm - cvsd
    lambda.min    <- mod.step2$lambda[which.min(cvm)]
    under.cvup1se <- which(cvm <= cvup[which.min(cvm)] * 1)
    min.nzero1se  <- under.cvup1se[which.min(mod.step2$nzero[under.cvup1se])]
    lambda.1se    <- mod.step2$lambda[min.nzero1se]
    
    naive.cv <- list()
    add.to.naive.cv <- c("cvm", "cvsd", "cvup", "cvlo", "lambda.min", 
                         "lambda.1se")
    for (add in add.to.naive.cv) {
      eval(parse(text = paste0("naive.cv$", add, " <- ", add)))}
    
    
    ada.lasso.naive.cv <- coef(mod.step2, s = naive.cv$lambda.min)
    
  }
  
  return(list(ada.lasso.naive.cv = ada.lasso.naive.cv,
              ada.lasso = mod.step2,
              naive.cv = naive.cv)
  )
}



