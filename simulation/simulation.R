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
## Required Packages: glmnet, Matrix, parallel, here, MASS, dplyr, ggplot2,
##                    openxlsx
## -----------------------------------------------------------------------------


set.seed(1)

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

source("../functions/ada.lasso.nested.cv.R")
source("../functions/ada.lasso.naive.cv.R")
source("../functions/evaluation.R")


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

if ( ! file.exists("intermediate_results/")) {
  dir.create("intermediate_results")
}
saveRDS(res.table, file = "intermediate_results/res_table.rds")






###################
###   Tables    ###
###################

### Table 1
# Averages of the computational times (in minutes)

tab.table <- res.table

tab.table$Runtime <- as.numeric(as.character(tab.table$Runtime))

main.tab.1 <- tab.table %>% group_by(Method) %>% summarise(Av = mean(Runtime)) 


# save Tables

if ( ! file.exists("tables/") ) {
  dir.create("tables")
}

write.xlsx(data.frame(main.tab.1), file = "tables/main_tab_1.xlsx",
           sheetName = "main_tab_1", append = FALSE)


###################
###   Figures   ###
###################


# formatting the results table

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



### Main Figure 2

main.fig.2 <-  
  ggplot(data = fig.table %>% filter(rho == 0.3,
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



### Supp Figure 1

supp.fig.1 <-  
  ggplot(data = fig.table %>% filter(rho == 0.6,
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



### Supp Figure 2

supp.fig.2 <-  
  ggplot(data = fig.table %>% filter(rho == 0.6,
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



### Supp Figure 3

supp.fig.3 <-  
  ggplot(data = fig.table %>% filter(rho == 0.6,
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



### Supp Figure 4

supp.fig.4 <-  
  ggplot(data = fig.table %>% filter(rho == 0.6,
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



### Supp Figure 5

supp.fig.5 <-  
  ggplot(data = fig.table %>% filter(rho == 0.3,
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



### Supp Figure 6

supp.fig.6 <-  
  ggplot(data = fig.table %>% filter(rho == 0.3,
                                     n == 500,
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


# save Figures

if ( ! file.exists("figures/") ) {
  dir.create("figures")
}

ggsave(filename = "figures/main_fig_2.pdf", plot = main.fig.2, 
       width = 8.5, height = 6, dpi = 600 , limitsize = FALSE, units = "in")




ggsave(filename = "figures/supp_fig_1.pdf", plot = supp.fig.1, 
       width = 8.5, height = 9, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "figures/supp_fig_2.pdf", plot = supp.fig.2, 
       width = 8.5, height = 9, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "figures/supp_fig_3.pdf", plot = supp.fig.3, 
       width = 8.5, height = 9, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "figures/supp_fig_4.pdf", plot = supp.fig.4, 
       width = 8.5, height = 9, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "figures/supp_fig_5.pdf", plot = supp.fig.5, 
       width = 8.5, height = 6, dpi = 600 , limitsize = FALSE, units = "in")

ggsave(filename = "figures/supp_fig_6.pdf", plot = supp.fig.6, 
       width = 8.5, height = 6, dpi = 600 , limitsize = FALSE, units = "in")




