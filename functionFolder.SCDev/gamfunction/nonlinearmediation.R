library(mgcv)
library(psych)
library(tidyverse)
library(mediation)
library(Formula)
library(stats)


#### FIT GAM SMOOTH ####
##Function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
nonlinearmediation <- function(mediator, dataname, cognition_var, smooth_var, covariates,knots, set_fx = FALSE){
  
  #Fit the gam
  gam.data <- get(dataname)
  NonNANIndex <- which(!is.na(gam.data[ ,cognition_var]))
  gam.data <- gam.data[NonNANIndex,]
  tmp<-gam.data[,mediator]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp), tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  parcel <- as.character(mediator)
  gam.data$tmp<-gam.data[ ,mediator]
  #fit models
  #X=smooth_var, M=mediator, Y=cognition_var
  modelformula_age <- as.formula(sprintf("%s ~s(%s, k = %s, fx = %s) + %s","tmp", smooth_var, knots, set_fx, covariates))
  med.fit <- gam(modelformula_age,method="REML", data = gam.data)
  modelformula_cog <- as.formula(sprintf("%s ~ %s +s(%s, k = %s, fx = %s) + %s",cognition_var, "tmp",smooth_var, knots, set_fx, covariates))
  out.fit <- gam(modelformula_cog,method="REML", data = gam.data)
  med.out<-mediate(med.fit, out.fit, treat = smooth_var, mediator = "tmp",
                   sims=1000, boot=TRUE)
  medparameter<-summary(med.out)
  indirect.Beta <- medparameter$d.avg
  indirect.lwth95CI <- medparameter$d.avg.ci[1]
  indirect.upth95CI <- medparameter$d.avg.ci[2]
  indirect.P <- medparameter$d.avg.p
  direct.Beta <- medparameter$z.avg
  direct.lwth95CI <- medparameter$z.avg.ci[1]
  direct.upth95CI <- medparameter$z.avg.ci[2]
  direct.P <- medparameter$z.avg.p
  total.Beta <- medparameter$tau.coef
  total.lwth95CI <- medparameter$tau.ci[1]
  total.upth95CI <- medparameter$tau.ci[2]
  total.P <- medparameter$tau.p
  prop.Med <- medparameter$n.avg
  prop.lwth95CI <- medparameter$n.avg.ci[1]
  prop.upth95CI <- medparameter$n.avg.ci[2]
  prop.P <- medparameter$n.avg.p
  
  stats.results <- cbind(parcel, cognition_var, indirect.Beta, indirect.lwth95CI, indirect.upth95CI,
                         indirect.P, direct.Beta, direct.lwth95CI, direct.upth95CI, direct.P, total.Beta,
                         total.lwth95CI, total.upth95CI, total.P, prop.Med, prop.lwth95CI, prop.upth95CI,
                         prop.P)
  return(stats.results)
}
  
  
  
  
  
  
  