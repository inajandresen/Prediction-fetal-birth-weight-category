
library(tidyverse)
library(Biobase)
library(splines)
library(mgcv)
library(parallel)

rm(list=ls())

load("../Data/winslgRFU_exprset_STORK_SGA_LGA.RData")
ano <- pData(winslgRFU_exprset_STORK_LGA_SGA)
exprs <- exprs(winslgRFU_exprset_STORK_LGA_SGA)
prot <- fData(winslgRFU_exprset_STORK_LGA_SGA)

ans=rownames(exprs)
for( i in 1:nrow(exprs)){
  an=ans[i]
  ano$Y=exprs[an,]
  #plot(Y~GAWeeks,data=ano[ano$Bwkat=="AGA",],pch=19,cex=0.6)
  ano$group=factor(ano$Bwkat)
  
  anoc<- ano[ano$group=="AGA",]
  anoc$ID<-factor(anoc$ID)

  fit3 <- gam(Y~s(GAWeeks,k=3)+s(ID,bs = 're'),
                          data=anoc,
                          method="ML")
  
  
  fit4 <- gam(Y~s(GAWeeks,k=4)+s(ID,bs = 're'),
                          data=anoc,
                          method="ML")
  
  fit5 <- gam(Y~s(GAWeeks,k=5)+s(ID,bs = 're'),
                          data=anoc,
                          method="ML")
  
  fit6 <- gam(Y~s(GAWeeks,k=6)+s(ID,bs = 're'),
                          data=anoc,
                          method="ML")
  
  fit=list(fit3,fit4,fit5,fit6)[[which.max(c(summary(fit3)$r.sq,
                                             summary(fit4)$r.sq,
                                             summary(fit5)$r.sq, 
                                             summary(fit6)$r.sq))]]
  
  ano1<-ano
  ano1$ID<-NULL
  ano$Ypred=predict.gam(fit,ano1,exclude="s(ID)", newdata.guaranteed=TRUE)
  ano$Yd=ano$Y-ano$Ypred
  #plot(Yd~GAWeeks,data=ano[ano$Bwkat=="AGA",],pch=19,cex=0.6)
  exprs[an,rownames(ano)] = ano$Yd
}


MoM_exprset_STORK_LGA_SGA <- ExpressionSet(exprs, phenoData = AnnotatedDataFrame(ano), featureData = AnnotatedDataFrame(prot))
save(MoM_exprset_STORK_LGA_SGA, file = "../Data/MoM_exprset_STORK_SGA_LGA.RData")
save(fit3, fit4, fit4, fit6, fit, file = "../Data/fit_MoM.RData")

