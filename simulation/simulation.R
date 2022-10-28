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
## Required Packages: glmnet, Matrix, parallel, MASS, dplyr, ggplot2, openxlsx
## -----------------------------------------------------------------------------


##########################
###   Read functions   ###
##########################

library(glmnet)
library(Matrix)
library(parallel)
library(MASS)
library(dplyr)
library(ggplot2)
library(openxlsx)


#### load instrumental functions
source("../functions/ada.lasso.nested.cv.R")
source("../functions/ada.lasso.naive.cv.R")
source("../functions/evaluation.R")


#### load auxialiary data for the computation of the correlation matrix inspired
#### by a real data example
load("for_simul_insp_real_data/order_metabolites.RData")
load("for_simul_insp_real_data/correlation_metabolites_all.RData")


####################
###   Settings   ###
####################

ncores     <- 25 #25 number of cores (used for the parallel computation below)
sim.seq    <- 1:100
signal.seq <- c(0.25, 0.5, 1, 1.5) / 4
eps.seq    <- c(0)
ntest      <- 10000


###############
###   Run   ###
###############

res.table=NULL
for (simtype in  c("Random", "RealInsp")) {
  if (simtype == "RealInsp") {
    rho.seq <- c(0.3)
  }# artificial, rho is not used in the Design inspired by the real data example
  if (simtype == "Random") {
    rho.seq <- c(0.3, 0.6)
  }
  for (rho in rho.seq) {
    if (rho == 0.3) {
      p.seq      <- c(100, 1000)
      p0.seq     <- c(10)
      ntrain.seq <- c(500)
    }
    if (rho == 0.6) {
      p.seq      <- c(100, 1000, 10000)
      p0.seq     <- c(5, 10, 50)
      ntrain.seq <- c(100, 500)
    }
    for(n in 1:length(ntrain.seq)) {
      ntot <- ntrain.seq[n] + ntest
      for(p in p.seq) {
        for(p0 in p0.seq) {
          for(signal in signal.seq) {
            for(eps in eps.seq) {
              
              funparall <- function (sim) {
                
                set.seed(sim)
                
                ev.table <- matrix(ncol = 6,nrow = 0)
                colnames(ev.table)=c("sACC","Precision","Recall","Pred Error",
                                     "Runtime", "Method")
                
                # data generation
                if (simtype == "RealInsp") {
                  ind_to_keep <- order.metabolites[1:p]
                  x <- mvrnorm(n = ntot, mu = rep(0, p),
                               Sigma = corr.mat[ind_to_keep, ind_to_keep])
                }
                if (simtype == "Random"){
                  x <- lapply(1:ntot, 
                              FUN = function (i) {
                                arima.sim(p, model = list(ar = rho))
                              })
                  x <- matrix(unlist(x), ncol = p, byrow = TRUE) 
                }
                
                beta <- rep(0, p)
                beta[sample(1:p, p0)] <- (2 * rbinom(n = p0, size = 1,
                                                     prob = 0.5) - 1) * signal
                
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
                
                opt.lasso <- coef(lasso, s = lasso$lambda.min)
                
                
                ada.lw.nai.start <- Sys.time()
                ada.lw.nai <- ada.lasso.naive.cv(x = xtrain, y = ytrain, 
                                                 alpha.weights = 1,
                                                 eps = eps, nfolds = 10)
                ada.lw.nai.end <- Sys.time()
                
                if (is.null(ada.lw.nai$ada.lasso.naive.cv)) {
                  opt.ada.lw.nai <- opt.lasso
                } else {
                  opt.ada.lw.nai <- ada.lw.nai$ada.lasso.naive.cv
                }
                
                
                ada.lw.nest.start <- Sys.time()
                ada.lw.nest <- ada.lasso.nested.cv(x = xtrain, y = ytrain, 
                                                   alpha.weights = 1,
                                                   eps = eps, nfolds = 10)
                ada.lw.nest.end <- Sys.time()
                
                if (is.null(ada.lw.nest$ada.lasso.nested.cv)) {
                  opt.ada.lw.nest <- opt.lasso 
                } else {
                  opt.ada.lw.nest <- ada.lw.nest$ada.lasso.nested.cv
                }
                
                if (rho == 0.3 & simtype == "Random") {
                  # adaptive lasso with ridge weights run only if 
                  # rho = 0.3 & simtype = "Random"
                  ada.rw.nai.start <- Sys.time()
                  ada.rw.nai <- ada.lasso.naive.cv(x = xtrain, y = ytrain, 
                                                   alpha.weights = 0,
                                                   eps = eps, nfolds = 10)
                  ada.rw.nai.end <- Sys.time()
                  
                  opt.ada.rw.nai <- ada.rw.nai$ada.lasso.naive.cv
                  
                  
                  ada.rw.nest.start <- Sys.time()
                  ada.rw.nest <- ada.lasso.nested.cv(x = xtrain, y = ytrain, 
                                                     alpha.weights = 0,
                                                     eps = eps, nfolds = 10)
                  ada.rw.nest.end <- Sys.time()
                  
                  opt.ada.rw.nest <- ada.rw.nest$ada.lasso.nested.cv
                } else {
                  ada.rw.nest <- ada.rw.nai <- NULL
                }
                
                if (p < ntrain.seq[n] & rho == 0.3 & simtype == "Random") {
                  # adaptive lasso with ols weigths run only if
                  # rho=0.3, p < n & simtype = "Random"
                  ada.ow.nai.start <- Sys.time()
                  ada.ow.nai <- ada.lasso.naive.cv(x = xtrain, y = ytrain, 
                                                   alpha.weights = 99,
                                                   eps = eps, nfolds = 10)
                  ada.ow.nai.end <- Sys.time()
                  
                  opt.ada.ow.nai <- ada.ow.nai$ada.lasso.naive.cv
                  
                  
                  ada.ow.nest.start <- Sys.time()
                  ada.ow.nest <- ada.lasso.nested.cv(x = xtrain, y = ytrain, 
                                                     alpha.weights = 99,
                                                     eps = eps, nfolds = 10)
                  ada.ow.nest.end <- Sys.time()
                  
                  opt.ada.ow.nest <- ada.ow.nest$ada.lasso.nested.cv
                } else {
                  ada.ow.nest <- ada.ow.nai <- NULL
                }
                
                ## save results for Fig 1
                
                if (simtype =="Random" & rho == 0.3 & ntrain.seq[n]  == 500 &
                    p == 1000 & signal == 1 / 4 & sim == 10) {
                  res.fig1 <<- list(data = list(xtrain = xtrain, xtest = xtest,
                                                ytrain = ytrain, ytest = ytest,
                                                beta = beta),
                                    res = list(lasso = lasso,
                                               ada.lw.nai = ada.lw.nai,
                                               ada.lw.nest = ada.lw.nest,
                                               ada.rw.nai = ada.rw.nai,
                                               ada.rw.nest = ada.rw.nest
                                    ))
                  if ( ! file.exists("intermediate_results/")) {
                    dir.create("intermediate_results")}
                  save(res.fig1, file = "intermediate_results/res_fig1.RData")
                }
                
                # Evaluation
                
                ev.table <- rbind(ev.table,
                                  c(evaluation(
                                    beta.star = beta,
                                    beta.hat = opt.lasso[-1],
                                    intercept.hat = opt.lasso[1], 
                                    y = ytest, x = xtest), 
                                    difftime(lasso.end, lasso.start,
                                             units = 'mins'), 
                                    "lasso"),
                                  c(evaluation(
                                    beta.star = beta, 
                                    beta.hat = opt.ada.lw.nai[-1], 
                                    intercept.hat = opt.ada.lw.nai[1],
                                    y = ytest, x = xtest), 
                                    difftime(ada.lw.nai.end, ada.lw.nai.start,
                                             units = 'mins'), 
                                    "adalasso.lw.nai"),
                                  c(evaluation(
                                    beta.star = beta, 
                                    beta.hat = opt.ada.lw.nest[-1],
                                    intercept.hat = opt.ada.lw.nest[1], 
                                    y = ytest, x = xtest),
                                    difftime(ada.lw.nest.end, ada.lw.nest.start,
                                             units = 'mins'),
                                    "adalasso.lw.nest"))
                
                
                if (rho == 0.3 & simtype == "Random") {
                  ev.table <- rbind(ev.table,            
                                    c(evaluation(
                                      beta.star = beta, 
                                      beta.hat = opt.ada.rw.nai[-1],
                                      intercept.hat = opt.ada.rw.nai[1], 
                                      y = ytest,x = xtest),
                                      difftime(
                                        ada.rw.nai.end, ada.rw.nai.start,
                                        units = 'mins'),
                                      "adalasso.rw.nai"),
                                    c(evaluation(
                                      beta.star = beta,
                                      beta.hat = opt.ada.rw.nest[-1], 
                                      intercept.hat = opt.ada.rw.nest[1], 
                                      y = ytest, x = xtest),
                                      difftime(
                                        ada.rw.nest.end, ada.rw.nest.start,
                                        units = 'mins'),
                                      "adalasso.rw.nest"))
                  
                  if (p < ntrain.seq[n]) {
                    ev.table <- rbind(ev.table,
                                      c(evaluation(
                                        beta.star = beta, 
                                        beta.hat = opt.ada.ow.nai[-1], 
                                        intercept.hat = opt.ada.ow.nai[1], 
                                        y = ytest, x = xtest),
                                        difftime(
                                          ada.ow.nai.end, ada.ow.nai.start,
                                          units = 'mins'),
                                        "adalasso.ow.nai"),
                                      c(evaluation(
                                        beta.star = beta, 
                                        beta.hat = opt.ada.ow.nest[-1], 
                                        intercept.hat = opt.ada.ow.nest[1],
                                        y = ytest, x = xtest),
                                        difftime(
                                          ada.ow.nest.end, ada.ow.nest.start,
                                          units = 'mins'),
                                        "adalasso.ow.nest")
                    )
                  }
                }
                
                ev.table <- data.frame(ev.table, rho = rho, eps = eps, p = p, 
                                       p0 = p0, signal = signal,
                                       n = ntrain.seq[n], sim = sim,
                                       simtype = simtype)
                
                return(ev.table)
              }
              
              res <- mclapply(sim.seq, FUN = funparall, mc.cores = ncores)
              
              sim.table <- NULL
              for (i in 1:length(sim.seq)) {
                sim.table <- rbind(sim.table, res[[i]])
              }
              
              res.table <- rbind(res.table, sim.table)
            }
          }
        }
      }
    }  
  }
}
# save results

if(!file.exists("intermediate_results/")){dir.create("intermediate_results")}
saveRDS(res.table,file = "intermediate_results/res_table.rds")



###################
##  Load results ##
###################

# load("intermediate_results/res_fig1.RData")
# res.table <- readRDS(("intermediate_results/res_table.rds"))


###################
###   Tables    ###
###################

### Table 1
# Averages of the computational times (in minutes)

tab.table <- res.table
tab.table$Runtime <- as.numeric(as.character(tab.table$Runtime))

main.tab.1 <- tab.table %>% group_by(simtype, rho, n, p, Method) %>% 
  summarise(Av = mean(Runtime, na.rm=T)) 

# save Tables

if ( ! file.exists("../tables/") ) {
  dir.create("../tables")
}

write.xlsx(data.frame(main.tab.1), file = "../tables/main_tab_1.xlsx",
           sheetName = "main_tab_1", append = FALSE)


###################
###   Figures   ###
###################

# Preparing of the results that illustrate the defect of the naive K-fold
# cross-validation scheme (section 2.3)

res  <- res.fig1$res
data <- res.fig1$data 

# Recovery of the prediction error estimated via the naive and the nested
# cross-validation scheme and calculation of the prediction error estimated
# on an independent test sample

fig1.table <- NULL
for (meth in names(res))
{
  if (meth == 'lasso') {
    lam.seq    <- res[[meth]]$lambda
    cvm.seq    <- res[[meth]]$cvm
    beta.mat   <- res[[meth]]$glmnet.fit$beta
    itcpt.seq  <- res[[meth]]$glmnet.fit$a0
    lam.x      <- lam.seq / max(lam.seq)
    
    fig1.table <- rbind(fig1.table, cbind(lam.x, cvm.seq, "cv nai", meth))
  } else {
    lam.seq    <- res[[meth]]$ada.lasso$lambda
    cvm.seq    <- res[[meth]][[3]]$cvm
    beta.mat   <- res[[meth]]$ada.lasso$beta
    itcpt.seq  <- res[[meth]]$ada.lasso$a0
    lam.x      <- lam.seq / max(lam.seq)
    
    fig1.table <- rbind(fig1.table, cbind(lam.x, cvm.seq,
                                          paste0("cv ", substr(meth, 8, 10)),
                                          substr(meth, 1, 6)))
  }
  
  ytest.mat <- matrix(rep(data$ytest, length(lam.x)), ncol = length(lam.x))
  itcpt.mat <- matrix(rep(itcpt.seq, length(data$ytest)), ncol = length(lam.x),
                      byrow=T)
  
  pred.error.seq <- colMeans((ytest.mat -
                                (data$xtest %*% beta.mat + itcpt.mat)) ^ 2)
  
  fig1.table <- rbind(fig1.table, cbind(lam.x, pred.error.seq,
                                        "test", substr(meth,1,6)))
}

# formatting the results table

fig1.table <- data.frame(fig1.table)
colnames(fig1.table) <- c("Lambda.ratio", "Pred.error", "Type", "Method")

fig1.table$Lambda.ratio <- as.numeric(as.character(fig1.table$Lambda.ratio))
fig1.table$Pred.error <- as.numeric(as.character(fig1.table$Pred.error))
fig1.table$Method  <- factor(fig1.table$Method, 
                             levels = c("lasso", "ada.lw", "ada.rw"), 
                             labels = c( "lasso", "one-step-lasso", 
                                         "ridge-ada-lasso"))
fig1.table$Type <- factor(fig1.table$Type, 
                          levels = c("test", "cv nai", "cv nes"), 
                          labels = c("Independent test sample", 
                                     '"Naive" CV', "Nested CV"))

# Selection of the value of the tuning parameter selected by each method and
# its value of the estimated true prediction error  

fig1.lines <- fig1.table %>% group_by(Method, Type) %>%
  summarise(Lambda.min = Lambda.ratio[which.min(Pred.error)])

fig1.lines <- bind_cols(fig1.lines, Pred.error = 
                          sapply(1:nrow(fig1.lines), function (i) {
                            fig1.table$Pred.error[fig1.table$Type == 
                                                    "Independent test sample" & 
                                                    fig1.table$Method == 
                                                    fig1.lines$Method[i] &
                                                    fig1.table$Lambda.ratio == 
                                                    fig1.lines$Lambda.min[i]][1]
                          }))



# Preparing of the simulation study results 

fig.table <- rbind(data.frame(Value = res.table$sACC,
                              Criterion = "sACC",
                              res.table[, -c(1:5)]),
                   data.frame(Value = res.table$Precision,
                              Criterion = "Precision",
                              res.table[, -c(1:5)]),
                   data.frame(Value = res.table$Recall,
                              Criterion = "Recall",
                              res.table[, -c(1:5)]),
                   data.frame(Value = res.table$Pred.Error,
                              Criterion = "Pred Error",
                              res.table[, -c(1:5)]))

# formatting the results table

fig.table$Value <- as.numeric(as.character(fig.table$Value))
fig.table$Criterion <- factor(fig.table$Criterion, levels = 
                                c("sACC", "Precision", "Recall", "Pred Error"))

fig.table$Method <- factor(fig.table$Method,
                           levels = c("adalasso.ow.nai",
                                      "adalasso.ow.nest",
                                      "adalasso.rw.nai",
                                      "adalasso.rw.nest",
                                      "adalasso.lw.nai",
                                      "adalasso.lw.nest",
                                      "lasso"),
                           labels = c("ols-ada lasso naive CV",
                                      "ols-ada lasso nested CV",
                                      "ridge-ada lasso naive CV",
                                      "ridge-ada lasso nested CV",
                                      "one-step lasso naive CV", 
                                      "one-step lasso nested CV",
                                      "lasso CV"))

# functions to use in ggplot()

Inf_95 <- function (x) {
  return(mean(x) - 1.96 * sd(x) / sqrt(length(x)))
}
Sup_95 <- function (x) {
  return(mean(x) + 1.96 * sd(x) / sqrt(length(x)))
}

# theme to use in ggplot()

theme <- theme(
  axis.text = element_text(size = 7),
  plot.background = element_rect(fill = "white"),
  strip.text.x = element_text(size = 12, face = "bold", color = "black"),
  strip.text.y = element_text(size = 12, face = "bold", color = "black"),
  text = element_text(size = 12,color = "black"),
  legend.position = "bottom",
  legend.title = element_blank(),
  legend.key = element_rect(fill = "white"),
  legend.background = element_rect(fill = "white"),
  legend.text = element_text(size = 12),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = "grey70",size = 0.1),
  panel.background = element_rect(fill = "white",color = "grey70", size = 0.1),
  strip.background = element_rect(colour = "grey30", fill = "white",size = 0.1))

col <- c("brown2", "chartreuse2", "blue")


### Main Figure 1

main.fig.1 <- 
  ggplot(data = fig1.table,
         aes(x = Lambda.ratio, y = Pred.error, colour = Type)) +
  geom_line(size = 1) + 
  facet_wrap( ~ Method, ncol = 3, scales = "free") +
  scale_x_log10(labels = function(x) sprintf("%g", x)) +
  ylim(0.95 * min(fig1.table$Pred.error), 1.8) +
  geom_vline(data = fig1.lines, aes(xintercept = Lambda.min, colour = Type),
             linetype = 2, size = 0.8, alpha = 0.7, show.legend = TRUE) + 
  geom_hline(data = fig1.lines, aes(yintercept = Pred.error, colour = Type),
             linetype = 1, size = 0.8, alpha = 0.7, show.legend = TRUE)+
  scale_color_manual(values = col[c(3,1,2)]) + theme + 
  xlab(expression(paste(Lambda, "/", Lambda, "_max"))) + ylab(" ") +
  guides(color = guide_legend(override.aes = list(size = 1.5)))


### Main Figure 2

main.fig.2 <-  
  ggplot(data = fig.table %>% filter(simtype == "Random",
                                     rho == 0.3,
                                     n == 500,
                                     Method %in% c("one-step lasso naive CV",
                                                   "one-step lasso nested CV",
                                                   "lasso CV")),
         mapping = aes(signal, Value, color = Method, group = Method)) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, linetype = 1) +
  stat_summary(fun = "mean", geom = "point", size = 2) +
  stat_summary(fun = "Inf_95", geom = "line", alpha = 0.7, linetype = 2,
               size = 0.3) +
  stat_summary(fun = "Sup_95",geom = "line", alpha =  0.7, linetype = 2,
               size = 0.3) +
  facet_grid( Criterion ~ p, scales = "free",
              labeller = labeller(Criterion = label_value, p = label_both)) +
  scale_color_manual(values = col) +
  scale_y_continuous(limits = function (b) {
    return(c(ifelse(min(b) == 1, 0.9, min(b)), max(b)))
  }) +
  scale_x_continuous(breaks = signal.seq,
                     expand = expansion(mult = c(0.1, 0.1))) + theme + 
  xlab("Signal strength") + ylab(" ") + 
  guides(color = guide_legend(override.aes = list(size = 5)))


### Main Figure 3

main.fig.3 <-  
  ggplot(data = fig.table %>% filter(simtype == "RealInsp",
                                     Method %in% c("one-step lasso naive CV",
                                                   "one-step lasso nested CV",
                                                   "lasso CV")),
         mapping = aes(signal, Value, color = Method, group = Method)) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, linetype = 1) +
  stat_summary(fun = "mean", geom = "point", size = 2) +
  stat_summary(fun = "Inf_95", geom = "line", alpha = 0.7, linetype = 2,
               size = 0.3) +
  stat_summary(fun = "Sup_95",geom = "line", alpha =  0.7, linetype = 2,
               size = 0.3) +
  facet_grid( Criterion ~ p, scales = "free",
              labeller = labeller(Criterion = label_value, p = label_both)) +
  scale_color_manual(values = col) +
  scale_y_continuous(limits = function (b) {
    return(c(ifelse(min(b) == 1, 0.9, min(b)), max(b)))
  }) +
  scale_x_continuous(breaks = signal.seq,
                     expand = expansion(mult = c(0.1, 0.1))) + theme + 
  xlab("Signal strength") + ylab(" ") + 
  guides(color = guide_legend(override.aes = list(size = 5)))


### Supp Figure 1

supp.fig.1 <-  
  ggplot(data = fig.table %>% filter(simtype == "Random",
                                     rho == 0.3,
                                     n == 500,
                                     Method %in% c("ridge-ada lasso naive CV",
                                                   "ridge-ada lasso nested CV",
                                                   "lasso CV")),
         mapping = aes(signal, Value, color = Method, group = Method)) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, linetype = 1) +
  stat_summary(fun = "mean", geom = "point", size = 2) +
  stat_summary(fun = "Inf_95", geom = "line", alpha = 0.7, linetype = 2,
               size = 0.3) +
  stat_summary(fun = "Sup_95",geom = "line", alpha =  0.7, linetype = 2,
               size = 0.3) +
  facet_grid( Criterion ~ p, scales = "free",
              labeller = labeller(Criterion = label_value, p = label_both)) +
  scale_color_manual(values = col) +
  scale_y_continuous(limits = function (b) {
    return(c(ifelse(min(b) == 1, 0.9, min(b)), max(b)))
  }) +
  scale_x_continuous(breaks = signal.seq,
                     expand = expansion(mult = c(0.1, 0.1))) + theme + 
  xlab("Signal strength") + ylab(" ") + 
  guides(color = guide_legend(override.aes = list(size = 5)))


### Supp Figure 2

supp.fig.2 <-  
  ggplot(data = fig.table %>% filter(simtype == "Random",
                                     rho == 0.3,
                                     n == 500,
                                     p == 100,
                                     Method %in% c("ols-ada lasso naive CV",
                                                   "ols-ada lasso nested CV",
                                                   "lasso CV")),
         mapping = aes(signal, Value, color = Method, group = Method)) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, linetype = 1) +
  stat_summary(fun = "mean", geom = "point", size = 2) +
  stat_summary(fun = "Inf_95", geom = "line", alpha = 0.7, linetype = 2,
               size = 0.3) +
  stat_summary(fun = "Sup_95",geom = "line", alpha =  0.7, linetype = 2,
               size = 0.3) +
  #facet_grid( Criterion ~ p, scales = "free",
  #            labeller = labeller(Criterion = label_value, p = label_both)) +
  facet_wrap( ~ p + Criterion, ncol = 2, scales = "free",
              labeller = labeller(Criterion = label_value, p = label_both)) +
  scale_color_manual(values = col) +
  scale_y_continuous(limits = function (b) {
    return(c(ifelse(min(b) == 1, 0.9, min(b)), max(b)))
  }) +
  scale_x_continuous(breaks = signal.seq,
                     expand = expansion(mult = c(0.1, 0.1))) + theme + 
  xlab("Signal strength") + ylab(" ") + 
  guides(color = guide_legend(override.aes = list(size = 5)))


### Supp Figure 3

supp.fig.3 <-  
  ggplot(data = fig.table %>% filter(simtype == "Random",
                                     rho == 0.6,
                                     p %in% c(100, 1000, 10000),
                                     Criterion == "Pred Error",
                                     Method %in% c("one-step lasso naive CV",
                                                   "one-step lasso nested CV")),
         mapping = aes(signal, Value, color = Method, group = Method)) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, linetype = 1) +
  stat_summary(fun = "mean", geom = "point", size = 2) +
  stat_summary(fun = "Inf_95", geom = "line", alpha = 0.7, linetype = 2,
               size = 0.3) +
  stat_summary(fun = "Sup_95",geom = "line", alpha =  0.7, linetype = 2,
               size = 0.3) +
  facet_grid( n + p0 ~ p, scales = "free",
              labeller = labeller(n = label_both, p0 = label_both,
                                  p = label_both)) +
  scale_color_manual(values = col[1:2]) +
  scale_y_continuous(limits = function (b) {
    return(c(ifelse(min(b) == 1, 0.9, min(b)), max(b)))
  }) +
  scale_x_continuous(breaks = signal.seq,
                     expand = expansion(mult = c(0.1, 0.1))) + theme + 
  xlab("Signal strength") + ylab(" ") + 
  guides(color = guide_legend(override.aes = list(size = 5)))


### Supp Figure 4

supp.fig.4 <-  
  ggplot(data = fig.table %>% filter(simtype == "Random",
                                     rho == 0.6,
                                     p %in% c(100, 1000, 10000),
                                     Criterion == "sACC",
                                     Method %in% c("one-step lasso naive CV",
                                                   "one-step lasso nested CV")),
         mapping = aes(signal, Value, color = Method, group = Method)) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, linetype = 1) +
  stat_summary(fun = "mean", geom = "point", size = 2) +
  stat_summary(fun = "Inf_95", geom = "line", alpha = 0.7, linetype = 2,
               size = 0.3) +
  stat_summary(fun = "Sup_95",geom = "line", alpha =  0.7, linetype = 2,
               size = 0.3) +
  facet_grid( n + p0 ~ p, scales = "free",
              labeller = labeller(n = label_both, p0 = label_both,
                                  p = label_both)) +
  scale_color_manual(values = col[1:2]) +
  scale_y_continuous(limits = function (b) {
    return(c(ifelse(min(b) == 1, 0.9, min(b)), max(b)))
  }) +
  scale_x_continuous(breaks = signal.seq,
                     expand = expansion(mult = c(0.1, 0.1))) + theme + 
  xlab("Signal strength") + ylab(" ") + 
  guides(color = guide_legend(override.aes = list(size = 5)))


### Supp Figure 5

supp.fig.5 <-  
  ggplot(data = fig.table %>% filter(simtype == "Random",
                                     rho == 0.6,
                                     p %in% c(100, 1000, 10000),
                                     Criterion == "Precision",
                                     Method %in% c("one-step lasso naive CV",
                                                   "one-step lasso nested CV")),
         mapping = aes(signal, Value, color = Method, group = Method)) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, linetype = 1) +
  stat_summary(fun = "mean", geom = "point", size = 2) +
  stat_summary(fun = "Inf_95", geom = "line", alpha = 0.7, linetype = 2,
               size = 0.3) +
  stat_summary(fun = "Sup_95",geom = "line", alpha =  0.7, linetype = 2,
               size = 0.3) +
  facet_grid( n + p0 ~ p, scales = "free",
              labeller = labeller(n = label_both, p0 = label_both,
                                  p = label_both)) +
  scale_color_manual(values = col[1:2]) +
  scale_y_continuous(limits = function (b) {
    return(c(ifelse(min(b) == 1, 0.9, min(b)), max(b)))
  }) +
  scale_x_continuous(breaks = signal.seq,
                     expand = expansion(mult = c(0.1, 0.1))) + theme + 
  xlab("Signal strength") + ylab(" ") + 
  guides(color = guide_legend(override.aes = list(size = 5)))


### Supp Figure 6

supp.fig.6 <-  
  ggplot(data = fig.table %>% filter(simtype == "Random",
                                     rho == 0.6,
                                     p %in% c(100, 1000, 10000),
                                     Criterion == "Recall",
                                     Method %in% c("one-step lasso naive CV",
                                                   "one-step lasso nested CV")),
         mapping = aes(signal, Value, color = Method, group = Method)) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, linetype = 1) +
  stat_summary(fun = "mean", geom = "point", size = 2) +
  stat_summary(fun = "Inf_95", geom = "line", alpha = 0.7, linetype = 2,
               size = 0.3) +
  stat_summary(fun = "Sup_95",geom = "line", alpha =  0.7, linetype = 2,
               size = 0.3) +
  facet_grid( n + p0 ~ p, scales = "free",
              labeller = labeller(n = label_both, p0 = label_both,
                                  p = label_both)) +
  scale_color_manual(values = col[1:2]) +
  scale_y_continuous(limits = function (b) {
    return(c(ifelse(min(b) == 1, 0.9, min(b)), max(b)))
  }) +
  scale_x_continuous(breaks = signal.seq,
                     expand = expansion(mult = c(0.1, 0.1))) + theme + 
  xlab("Signal strength") + ylab(" ") + 
  guides(color = guide_legend(override.aes = list(size = 5)))


# save Figures

if ( ! file.exists("../figures/") ) {
  dir.create("../figures")
}

ggsave(filename = "../figures/main_fig_1.pdf", plot = main.fig.1, 
       width = 8.5, height = 4, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "../figures/main_fig_2.pdf", plot = main.fig.2, 
       width = 6.5, height = 6, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "../figures/main_fig_3.pdf", plot = main.fig.3, 
       width = 6.5, height = 6, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "../figures/supp_fig_1.pdf", plot = supp.fig.1, 
       width = 6.5, height = 6, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "../figures/supp_fig_2.pdf", plot = supp.fig.2, 
       width = 6.5, height = 6, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "../figures/supp_fig_3.pdf", plot = supp.fig.3, 
       width = 8.5, height = 9, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "../figures/supp_fig_4.pdf", plot = supp.fig.4, 
       width = 8.5, height = 9, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "../figures/supp_fig_5.pdf", plot = supp.fig.5, 
       width = 8.5, height = 9, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "../figures/supp_fig_6.pdf", plot = supp.fig.6, 
       width = 8.5, height = 9, dpi = 600 , limitsize = FALSE, units = "in")
