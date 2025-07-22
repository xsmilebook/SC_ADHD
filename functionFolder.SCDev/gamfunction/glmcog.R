library(tidyr)
library(mgcv)
library(psych)
library(tidyverse)
library(dplyr)


#### FIT GAM SMOOTH ####
##Function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
glm.fit.cognition <- function(region, dataname, cognition_var, smooth_var, covariates,corrmethod, stats_only = FALSE){
  
  #Fit the glm
  glm.data <- get(dataname)
  cognition<-glm.data[ ,cognition_var]
  outlierindx1<-which(cognition<mean(cognition,na.rm = T)-3*sd(cognition,na.rm = T) | cognition>mean(cognition,na.rm = T)+3*sd(cognition,na.rm = T))
  glm.data[outlierindx1 ,cognition_var]<-NA
  NonNANIndex <- which(!is.na(glm.data[ ,cognition_var]) & !is.na(glm.data[ ,region]))
  
  glm.data <- glm.data[NonNANIndex,]
  tmp<-glm.data[,region]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    glm.data<-glm.data[-outlierindx, ]
  }
  parcel <- as.character(region)
  modelformula <- as.formula(sprintf("%s ~ %s + %s + %s",region, cognition_var, smooth_var, covariates))
  modelformula.null<-as.formula(sprintf("%s ~ %s + %s",region, smooth_var, covariates))
  glm.model <- glm(modelformula, method="glm.fit", data = glm.data)
  glm.results <- summary(glm.model)
  glm.model.null <- glm(modelformula.null, method="glm.fit", data = glm.data)
  
  #glm statistics
  #F value for the smooth term and glm-based significance of the smooth term
  glm.smooth.t <- glm.results$coefficients[2,3]
  glm.smooth.pvalue <- glm.results$coefficients[2,4]
  
  #Full versus reduced model anova p-value
  anova.cov.pvalue <- anova(glm.model.null,glm.model,test='Chisq')$`Pr(>Chi)`[2]
  if(is.na(anova.cov.pvalue)){ #if residual deviance is exactly equal between full and reduced models and p=value = NA, set p = 1
    anova.cov.pvalue <- 1}
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((glm.model$y - glm.model$fitted.values)^2)
  sse.nullmodel <- sum((glm.model.null$y - glm.model.null$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  ### effect direction
  if(glm.smooth.t < 0){ #if the glm t-value for covariate of interest is less than 0, make the partialRsq negative
    partialRsq <- partialRsq*-1}
  
  #residual correlation
  varcorformula1 <- as.formula(sprintf("%s ~ %s + %s", region, smooth_var, covariates))
  varcorformula2 <- as.formula(sprintf("%s ~ %s + %s", cognition_var, smooth_var, covariates))
  res1<-residuals(glm(varcorformula1,method="glm.fit", data = glm.data))
  res2<-residuals(glm(varcorformula2,method="glm.fit", data = glm.data))
  PCorr_Test <- corr.test(res1, res2, method=corrmethod)
  correstimate<-as.numeric(PCorr_Test$r)
  corrp <- as.numeric(PCorr_Test$p)
  
  stats.results <- cbind(parcel, cognition_var, glm.smooth.t, glm.smooth.pvalue,anova.cov.pvalue, partialRsq, correstimate, corrp)
  data.results <- list()
  data.results[[1]] <- as.data.frame(stats.results)
  data.results[[2]] <- data.frame(SCres=res1, cogres=res2)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(data.results)
}





