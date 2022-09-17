## -----------------------------------------------------------------------------
## Filename: simulation.R
## -----------------------------------------------------------------------------
## Authors: Nadim Ballout, Lola Etievant, Vivian Viallon
## -----------------------------------------------------------------------------
## BiomJ article: "On the use of cross-validation for the calibration of
##                 the adaptive lasso"
## -----------------------------------------------------------------------------
## Sections in the article: main(section 2.3 and 4) + suppmat
## -----------------------------------------------------------------------------
## Required Packages: glmnet, Matrix, parallel, here, MASS
## -----------------------------------------------------------------------------


set.seed(1)

##########################
###   Read functions   ###
##########################

library(glmnet)
library(Matrix)
library(parallel)
library(MASS)

source("functions/ada.lasso.nested.cv.R")
source("functions/ada.lasso.naive.cv.R")
source("functions/evaluation.R")


####################
###   Settings   ###
####################

ncores     <- 1 # number of cores (used for the parallel computation below)
sim.seq    <- 1:100
signal.seq <- c(0.25, 0.5, 1, 1.5) / 4
eps.seq    <- c(0)
p.seq      <- c(100, 500, 1000, 10000)
p0.seq     <- c(5, 10, 50)
ntrain.seq <- c(100, 500)
ntest      <- 10000
rho.seq    <- c(0.3, 0.6)

###############
###   Run   ###
###############

res.table=NULL

for (rho in rho.seq) {
  for(n in 1:length(ntrain.seq)){
    ntot=ntrain.seq[n]+ntest
    for(p in p.seq){
      for(p0 in p0.seq){
        for(signal in signal.seq){
          for(eps in eps.seq){
            
            funparall <- function (sim) {
              
              ev.table <- matrix(ncol = 6,nrow = 0)
              colnames(ev.table)=c("sACC","Precision","Recall","Pred Error",
                                   "Runtime","Method")
              
              # data generation
              
              x <- lapply(1:ntot, 
                          FUN = function(i){arima.sim(p, model =
                                                        list(ar = rho))})
              
              x <- matrix(unlist(x), ncol = p, byrow = TRUE) 
              
              beta <- rep(0, p)
              beta[sample(1:p, p0)] <- (2 * rbinom(n = p0, size = 1, prob = 0.5)
                                        - 1) * signal
              
              y <- x %*% beta + rnorm(ntot, 0, 1)
              
              xtrain <- x[1:ntrain.seq[n], ]
              ytrain <- y[1:ntrain.seq[n]]
              
              xtest <- x[(ntrain.seq[n] + 1):ntot, ]
              ytest <- y[(ntrain.seq[n] + 1):ntot]
              
              # coefficients estimation
              
              lasso.start <- Sys.time()
              lasso <- cv.glmnet(x = xtrain, y = ytrain, family = "gaussian",  
                                 foldid = rep_len(1:10, length(ytrain)),
                                 alpha = 1, standardize = TRUE)
              lasso.end <- Sys.time()
              
              opt.lasso <- coef(lasso,s=lasso$lambda.min)
              
              
              ada.lw.nai.start <- Sys.time()
              ada.lw.nai <- ada.lasso.naive.cv(x = xtrain, y = ytrain, 
                                               alpha.weights = 1, eps = eps, 
                                               nfolds = 10)
              ada.lw.nai.end <- Sys.time()
              
              if (is.null(ada.lw.nai$ada.lasso.naive.cv)) {
                opt.ada.lw.nai <- opt.lasso
              } else {
                opt.ada.lw.nai <- ada.lw.nai$ada.lasso.naive.cv
              }
              
              
              ada.lw.nest.start <- Sys.time()
              ada.lw.nest <- ada.lasso.nested.cv(x = xtrain, y = ytrain, 
                                                 alpha.weights = 1, eps = eps, 
                                                 nfolds = 10)
              ada.lw.nest.end <- Sys.time()
              
              if (is.null(ada.lw.nest$ada.lasso.nested.cv)) {
                opt.ada.lw.nest <- opt.lasso 
              } else {
                opt.ada.lw.nest <- ada.lw.nest$ada.lasso.nested.cv
              }
              
              
              ada.rw.nai.start <- Sys.time()
              ada.rw.nai <- ada.lasso.naive.cv(x = xtrain, y = ytrain, 
                                               alpha.weights = 0, eps = eps, 
                                               nfolds = 10)
              ada.rw.nai.end <- Sys.time()
              
              opt.ada.rw.nai <- ada.rw.nai$ada.lasso.naive.cv
              
              
              ada.rw.nest.start <- Sys.time()
              ada.rw.nest <- ada.lasso.nested.cv(x = xtrain, y = ytrain, 
                                                 alpha.weights = 0, eps = eps, 
                                                 nfolds = 10)
              ada.rw.nest.end <- Sys.time()
              
              opt.ada.rw.nest <- ada.rw.nest$ada.lasso.nested.cv
              
              
              if (p < ntrain.seq[n]) {
                
                ada.ow.nai.start <- Sys.time()
                ada.ow.nai <- ada.lasso.naive.cv(x = xtrain, y = ytrain, 
                                                 alpha.weights = 99, eps = eps, 
                                                 nfolds = 10)
                ada.ow.nai.end <- Sys.time()
                
                opt.ada.ow.nai <- ada.ow.nai$ada.lasso.naive.cv
                
                
                ada.ow.nest.start <- Sys.time()
                ada.ow.nest <- ada.lasso.nested.cv(x = xtrain, y = ytrain, 
                                                   alpha.weights = 99, eps = eps, 
                                                   nfolds = 10)
                ada.ow.nest.end <- Sys.time()
                
                opt.ada.ow.nest <- ada.ow.nest$ada.lasso.nested.cv
                
                
              } else {
                
                ada.ow.nest <- ada.ow.nai <- NULL
                
              }
              
              # Evaluation
              
              ev.table <- rbind(ev.table,
                                c(evaluation(beta.star = beta,
                                             beta.hat = opt.lasso[-1],
                                             intercept.hat = opt.lasso[1], 
                                             y = ytest, x = xtest), 
                                  lasso.end - lasso.start,
                                  "lasso"),
                                c(evaluation(beta.star = beta, 
                                             beta.hat = opt.ada.lw.nai[-1], 
                                             intercept.hat = opt.ada.lw.nai[1],
                                             y = ytest, x = xtest), 
                                  ada.lw.nai.end - ada.lw.nai.start, 
                                  "adalasso.lw.nai"),
                                c(evaluation(beta.star = beta, 
                                             beta.hat = opt.ada.lw.nest[-1],
                                             intercept.hat = opt.ada.lw.nest[1], 
                                             y = ytest, x = xtest),
                                  ada.lw.nest.end - ada.lw.nest.start,
                                  "adalasso.lw.nest"),
                                c(evaluation(beta.star = beta, 
                                             beta.hat = opt.ada.rw.nai[-1],
                                             intercept.hat = opt.ada.rw.nai[1], 
                                             y = ytest,x = xtest),
                                  ada.rw.nai.end - ada.rw.nai.start,
                                  "adalasso.rw.nai"),
                                c(evaluation(beta.star = beta,
                                             beta.hat = opt.ada.rw.nest[-1], 
                                             intercept.hat = opt.ada.rw.nest[1], 
                                             y = ytest, x = xtest),
                                  ada.rw.nest.end - ada.rw.nest.start,
                                  "adalasso.rw.nest"))
              
              if (p < ntrain.seq[n]) {
                ev.table <- rbind(ev.table,
                                  c(evaluation(beta.star = beta, 
                                               beta.hat = opt.ada.ow.nai[-1], 
                                               intercept.hat = opt.ada.ow.nai[1], 
                                               y = ytest, x = xtest),
                                    ada.ow.nai.end - ada.ow.nai.start,
                                    "adalasso.ow.nai"),
                                  c(evaluation(beta.star = beta, 
                                               beta.hat = opt.ada.ow.nest[-1], 
                                               intercept.hat = opt.ada.ow.nest[1],
                                               y = ytest, x = xtest),
                                    ada.ow.nest.end - ada.ow.nest.start,
                                    "adalasso.ow.nest")
                )  
              }
              
              ev.table <- data.frame(ev.table, rho = rho, eps = eps, p = p, 
                                     p0 = p0, signal = signal,
                                     n = ntrain.seq[n], sim = sim)
              
              return(ev.table)
            }
            
            res <- mclapply(sim.seq, FUN = funparall, mc.cores = ncores)
            
            sim.table <- NULL
            for (i in 1:length(sim.seq)) {
              sim.table <- rbind(sim.table, res[[i]])}
            
            res.table <- rbind(res.table, sim.table)
          }
        }
      }
    }
  }  
}

# save results

save(res.table,file = "simulation/intermediate_results/res.table.Rdata")


###################
###   Figures   ###
###################

