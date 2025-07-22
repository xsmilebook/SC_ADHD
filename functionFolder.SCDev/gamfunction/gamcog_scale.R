library(tidyr)
library(mgcv)
library(psych)
library(tidyverse)
library(dplyr)


#### FIT GAM SMOOTH ####
# for cognition score which is adjusted for age
##Function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
gam.fit.cognition <- function(region, dataname, cognition_var, smooth_var, covariates, knots,corrmethod, set_fx = FALSE, stats_only = FALSE){
  
  #Fit the gam
  gam.data <- get(dataname)
  cognition<-gam.data[ ,cognition_var]
  outlierindx1<-which(cognition<mean(cognition,na.rm = T)-3*sd(cognition,na.rm = T) | 
                      cognition>mean(cognition,na.rm = T)+3*sd(cognition,na.rm = T))
  gam.data[outlierindx1 ,cognition_var]<-NA
  NonNANIndex <- which(!is.na(gam.data[ ,cognition_var]))
  
  gam.data <- gam.data[NonNANIndex,]
  tmp<-gam.data[,region]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  parcel <- as.character(region)
  nnum<-nrow(gam.data)
  
  parcel <- as.character(region)
  modelformula <- as.formula(sprintf("%s ~ %s + %s",region, cognition_var, covariates))
  modelformula.null<-as.formula(sprintf("%s ~ %s",region, covariates))
  gam.model <- gam(modelformula, method="REML", data = gam.data)
  gam.results <- summary(gam.model)
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data)
  
  #GAM statistics
  #F value for the smooth term and GAM-based significance of the smooth term
  gam.smooth.t <- gam.results$p.table[2,3]
  gam.smooth.pvalue <- gam.results$p.table[2,4]
  
  #Full versus reduced model anova p-value
  anova.cov.pvalue <- anova.gam(gam.model.null,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  if(is.na(anova.cov.pvalue)){ #if residual deviance is exactly equal between full and reduced models and p=value = NA, set p = 1
    anova.cov.pvalue <- 1}
  
  #residual correlation
  varcorformula1 <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  varcorformula2 <- as.formula(sprintf("%s ~ %s", cognition_var, covariates))
  res1<-residuals(gam(varcorformula1,method="REML", data = gam.data))
  res2<-residuals(gam(varcorformula2,method="REML", data = gam.data))
  PCorr_Test <- corr.test(res1, res2, method=corrmethod)
  correstimate<-PCorr_Test$r
  corrp <- PCorr_Test$p
  
  stats.results <- cbind(parcel,nnum, cognition_var, gam.smooth.t, gam.smooth.pvalue,anova.cov.pvalue, correstimate, corrp)
  data.results <- list()
  data.results[[1]] <- as.data.frame(stats.results)
  data.results[[2]] <- data.frame(SCres=res1, cogres=res2)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(data.results)
}





