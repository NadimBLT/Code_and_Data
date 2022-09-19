## -----------------------------------------------------------------------------
## Title: R-function evaluation()
## -----------------------------------------------------------------------------
## Authors: Nadim Ballout, Lola Etievant, Vivian Viallon
## -----------------------------------------------------------------------------
## Description:  This function calculates and returns the (signed) support 
##               accuracy, the precision, the recall and the prediction error
## -----------------------------------------------------------------------------
## BiomJ article: "On the use of cross-validation for the calibration of
##                 the adaptive lasso"
## -----------------------------------------------------------------------------
## Section in the article: 4
## -----------------------------------------------------------------------------
## Required Package: -
## -----------------------------------------------------------------------------
## Usage: evaluation(
##                   x,
##                   y,
##                   beta.star,
##                   beta.hat,
##                   intercept.hat
##                   )
## -----------------------------------------------------------------------------
## Arguments:
##
##  x                input matrix, of dimension n x p, where n is the number of
##                   observations and p is the number of variables; each row is
##                   an observation vector.
##
##  y                quantitative response variable.
##
##  beta.star        vector, of length p, containing the true beta values
##
##  beta.hat         vector, of length p, containing the estimated beta values
##
##  intercept.hat    the estimated intercept value
## -----------------------------------------------------------------------------
## Value:            vector, of length 4, containing the value of the (signed) 
##                   support accuracy, the precision, the recall and the
##                   prediction error
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

evaluation <- function (beta.star, beta.hat, intercept.hat, x, y) {
  
  pred.error <- mean((y - (x %*% beta.hat + intercept.hat)) ^ 2)
  
  precision <- length(which(beta.hat != 0 & beta.star != 0)) / length(which(
    beta.hat  != 0))
  
  recall    <- length(which(beta.hat != 0 & beta.star != 0)) / length(which(
    beta.star != 0))
  
  beta.star[which(beta.star<0)] <- -1
  beta.star[which(beta.star>0)] <-  1
  beta.hat[which(beta.hat<0)]   <- -1
  beta.hat[which(beta.hat>0)]   <-  1
  
  sACC  <- length(which(beta.hat == beta.star)) / length(beta.star) -
    length(which(beta.hat * beta.star == -1)) / length(beta.star)
  
  vec <- c(sACC, precision, recall, pred.error)
  names(vec) <- c("sACC", "Precision", "Recall", "Pred Error")
  
  return(vec)
}



