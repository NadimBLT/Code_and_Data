rm(list=ls())
library(Matrix)
library(glmnet)
library(parallel)
library(here)
library(ggplot2)
library(tidyverse)

#path="C://Users/ballout/Desktop/CDESS/PAPIER-3/Codes/N/Codes/"

source("AdapGlmnet-NestedCv.R")

#settings

set.seed(135490)

vp      = c(100,500,1000)
vsignal = c(0.25,0.5,1,2)
vs0     = c(10, 50)
eps     = 1e-4
ntrain  = c(1000,2000)
ntest   = 10000


code.meth=rbind(c(FALSE,rep(TRUE,3)),c(NA,1,0,99))
colnames(code.meth)=c("Lasso","Onestep","Adlasso.ridge","Adlasso.mco")
rownames(code.meth)=c("adap","alpha.weights")

ntrain = 1000; p=1000; s0=10; signal = 0.5
res=list()

ntot=ntrain + ntest
x=matrix(rnorm(p*ntot),ncol = p)

beta=rep(0,p)
beta[sample(1:p,s0)]= (2*rbinom(n = s0,size = 1,prob = 0.5)-1)*signal

y=x%*%beta+rnorm(ntot,0,1)

xtrain=x[1:ntrain,]
ytrain=y[1:ntrain]

xtest=x[(ntrain+1):ntot,]
ytest=y[(ntrain+1):ntot]

if(p>=(ntrain/2)){valpha.weights=c(1,0)}else{valpha.weights=c(1,0,99)}


for (adap in c(FALSE,TRUE)) {
if (adap==TRUE){
  for (alpha.weights in valpha.weights) {
    
    mod=adap.glmnet(adap=adap,alpha.weights=alpha.weights,eps=eps,x=xtrain,y=ytrain,family="gaussian")#,alpha =1, intercept = TRUE)
    
    eval(parse(text=paste0("res$",names(which(code.meth[1,]==adap & code.meth[2,]==alpha.weights)),"=mod")))}
}else{
  
  mod=adap.glmnet(adap=adap,x=xtrain,y=ytrain,family="gaussian")#,alpha =1, intercept = TRUE)
  
  res$Lasso=mod}
}

RECAP_PREDS <- NULL
for (nam in names(res))
{
  lamtemp = res[[nam]]$lambda
  cvmtemp = res[[nam]]$cvm
  RECAP_PREDS <- rbind(RECAP_PREDS, cbind(lamtemp, cvmtemp, rep("CV", length(lamtemp)), rep(nam, length(lamtemp))))
  if (nam !='Lasso')
  {
    nestedcvm <- res[[nam]]$nested.cvm 
    RECAP_PREDS <- rbind(RECAP_PREDS, cbind(lamtemp, nestedcvm, rep("nestedCV", length(lamtemp)), rep(nam, length(lamtemp)) ))
  }
  
  BetaTemp <- res[[nam]]$glmnet.fit$beta
  IndTemp  <- res[[nam]]$glmnet.fit$a0
  
  MatYTest <- matrix(rep(ytest, length(lamtemp)), ncol=length(lamtemp))
  MatIntTe <- matrix(rep(IndTemp, length(ytest)), ncol=length(lamtemp), byrow=T) 
  PredErrTemp <- colMeans((MatYTest - (xtest%*%BetaTemp + MatIntTe ))^2)
  RECAP_PREDS <- rbind(RECAP_PREDS, cbind(lamtemp, PredErrTemp, rep("Test", length(lamtemp)), rep(nam, length(lamtemp))))
}

RECAP_PREDS <- data.frame(RECAP_PREDS)
colnames(RECAP_PREDS) <- c("Lambda", "PredError", "Type", "Method")
RECAP_PREDS$Lambda <- as.numeric(as.character(RECAP_PREDS$Lambda))
RECAP_PREDS$PredError <- as.numeric(as.character(RECAP_PREDS$PredError))
levels(RECAP_PREDS$Method) = c("Ridge-AdaLasso", "Lasso", "1-step-Laaasso") #"OLS-AdaLasso", 
RECAP_PREDS$Method <- factor(RECAP_PREDS$Method, levels= levels(RECAP_PREDS$Method)[c(2, 3, 1)]) #[c(3, 4, 1, 2)]


DataMins    <- RECAP_PREDS %>% group_by(Method, Type) %>% summarise( LambdaMin=Lambda[which.min(PredError)])


ggplot(data=RECAP_PREDS , aes(x=Lambda, y=PredError, colour=Type)) + geom_line() +  #%>% filter(Method != "Ridge-AdaLasso")
  facet_wrap(~Method, ncol=4, scales="free") + scale_x_log10() + #ylim(0.95*min(RECAP_PREDS$PredError), 2) + 
  theme(legend.position = "bottom", legend.title = element_blank(), axis.title.y = element_blank()) + 
  geom_vline(data=DataMins, aes(xintercept=LambdaMin, colour=Type),linetype="dotted")

ggsave(file="FIGs/ForFig1_new.png", width=6, height = 3.5, units="in", dpi=300)



