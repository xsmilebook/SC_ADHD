## This script is a function which fits linear mixed models with random intercepts, and compute the 
## individual slopes.

library(tidyverse)
library(lme4)

lmm.idv.randslope <- function(dep.var, dataname, indep.var, covariates = NA, subvar, mod_only = FALSE){
  
  #Fit the lm
  lm.data <- get(dataname)
  parcel <- as.character(dep.var)
  if (is.na(covariates)){
    modelformula <- as.formula(sprintf("%1$s ~ 1 + %2$s + (1 + %2$s | %3$s)", dep.var, indep.var, subvar))
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ 1 + %2$s + %3$s + (1 + %2$s | %4$s)", dep.var, indep.var, covariates, subvar))
  }
  
  lm.mod <- lmer(modelformula, data = lm.data, control = lmerControl(calc.derivs = FALSE),REML=TRUE)
  fixed=summary(lm.mod)$coefficient
  rand_int=ranef(lm.mod)[[subvar]]
  
  
  slopedf <- data.frame(subID =rownames(rand_int), slope = rand_int[[indep.var]]+fixed[2,1])
  names(slopedf) <- c("subID", paste0("slope.", dep.var))
  
  
  return(slopedf)
}





