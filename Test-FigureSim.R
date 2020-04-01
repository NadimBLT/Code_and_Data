load("TableALL-NewTest-5-simulations.Rdata")
View(TableALL)

TableALLP=rbind(as.matrix(cbind(TableALL[,-1],"PredError")),as.matrix(cbind(TableALL[,-2],"Accuracy")))
TableALLP=as.data.frame(TableALLP)
names(TableALLP)[1]="Value"
names(TableALLP)[9]="Criterion"
TableALLP$Approach=as.factor(paste0(TableALLP$Method,"-",TableALLP$Type))
TableALLP$Value=as.numeric(as.character(TableALLP$Value))
TableALLP$Approach=factor(TableALLP$Approach,levels = c("OLS.AdaLasso-CV","OLS.AdaLasso-nestedCV","One.step.Lasso-CV","One.step.Lasso-nestedCV","Ridge.AdaLasso-CV","Ridge.AdaLasso-nestedCV","Lasso-CV"))

Approach=c("OLS.AdaLasso-CV","OLS.AdaLasso-nestedCV","One.step.Lasso-CV","One.step.Lasso-nestedCV","Ridge.AdaLasso-CV","Ridge.AdaLasso-nestedCV","Lasso-CV")

Col=c("darkgoldenrod4","darkgoldenrod1","brown4","brown1","chartreuse4","chartreuse1","blue")

names(Col)=Approach


library(ggplot2)


ggplot(data = TableALLP,mapping = aes(signal,Value,color=Approach,group=Approach))+
  stat_summary(fun.y="mean",geom = "line",size=0.5,linetype = 1) +
  stat_summary(fun.y="mean",geom = "point",size=1.5) +
  facet_grid(n*Criterion ~  s0*p   ,scales = "free")+
  scale_color_manual(values = Col,name=" ")+
  # scale_y_continuous(name="Mesure") + 
  #scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
  guides(linetype = guide_legend(order = 1), color = guide_legend(order = 2))+
  theme( axis.text.x = element_text(hjust = 0.4),
         #axis.text.x=element_blank(),
         # axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="bottom",
         legend.key = element_rect(fill = "white"),
         legend.background = element_rect(fill = "white"),
         #legend.position = c(0.14, 0.80),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_line(colour = "grey95",size = 0.1),
         panel.background = element_rect(fill = "white",color = "grey50", size = 0.1),
         strip.background = element_rect(colour = "grey50", fill = "white",size = 0.1)) 


Inf_95 <- function(x){return(mean(x)-1.96*sd(x)/sqrt(length(x)))}
Sup_95 <- function(x){return(mean(x)+1.96*sd(x)/sqrt(length(x)))}



ggplot(data = TableALLP[-which(TableALLP$Method=="Ridge.AdaLasso"),],mapping = aes(signal,Value,color=Approach,group=Approach))+
  stat_summary(fun.y="mean",geom = "line",size=0.5,linetype = 1) +
  stat_summary(fun.y="mean",geom = "point",size=1.5) +
  stat_summary(fun.y="Inf_95",geom = "line",alpha=0.9,linetype=2,size=0.25) +
  stat_summary(fun.y="Sup_95",geom = "line",alpha=0.9,linetype=2,size=0.25) +
  facet_grid(n*Criterion ~  s0*p   ,scales = "free", labeller = "label_both")+
  scale_color_manual(values = Col[-c(5,6)],name=" ")+
  # scale_y_continuous(name="Mesure") + 
  #scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
  guides(linetype = guide_legend(order = 1), color = guide_legend(order = 2))+
  theme( axis.text.x = element_text(hjust = 0.4),
         #axis.text.x=element_blank(),
         # axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="bottom",
         legend.key = element_rect(fill = "white"),
         legend.background = element_rect(fill = "white"),
         #legend.position = c(0.14, 0.80),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_line(colour = "grey95",size = 0.1),
         panel.background = element_rect(fill = "white",color = "grey50", size = 0.1),
         strip.background = element_rect(colour = "grey50", fill = "white",size = 0.1)) 


