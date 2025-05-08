rm(list=ls())
library(tidyverse)
library(Biobase)
library(reshape2)
library(ggforce)
#library(pcaMethods)
library(scales)
library(gridExtra)
library(ggpubr)

#Make the function 
winsorize<-function(x)
{
  if(any(x> quantile(x,0.98)*2)){
    x[x>quantile(x,0.98)*2] <- quantile(x,0.98)*2
  }
  return(x)
  
}

#Outlier are defined as any data point above the 2*98th quantile (for that specific protein). 
#The winsorization moves outliers to the 2*98th quantile

load("../Data/RFU_exprset_STORK_SGA_LGA.RData")
exprs <- exprs(RFU_exprset_STORK_LGA_SGA)
samp <- pData(RFU_exprset_STORK_LGA_SGA)
prot <- fData(RFU_exprset_STORK_LGA_SGA)

#Are there any out lier?
outliers<- apply(exprs,1,function(x){any(x> quantile(x,0.98)*2)})
table(outliers)
#FALSE  TRUE 
#3096  1469 

exprs_winsorized <- t(apply(exprs,1,winsorize))
par(mfrow= c(2,1))
plot(density(exprs["CRYBB2.10000.28",]),main="Before")
plot(density(exprs_winsorized["CRYBB2.10000.28",]),main="After")

#Log2 tranform
exprs_log <- log2(exprs)
exprs_winsorized_log <- log2(exprs_winsorized)
par(mfrow= c(2,1))
plot(density(exprs_log["CRYBB2.10000.28",]),main="Before")
plot(density(exprs_winsorized_log["CRYBB2.10000.28",]),main="After")


#Make and save expressionset, winslgRFU

winslgRFU_exprset_STORK_LGA_SGA <- ExpressionSet(exprs_winsorized_log, phenoData = AnnotatedDataFrame(samp), featureData = AnnotatedDataFrame(prot))
save(winslgRFU_exprset_STORK_LGA_SGA, file = "../Data/winslgRFU_exprset_STORK_SGA_LGA.RData")

