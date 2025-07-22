library(lme4)
library(gamm4)
library(tidyr)
library(mgcv)
library(psych)
library(tidyverse)
library(dplyr)
library(lavaan)

set.seed(0925)
linearmediation_longitudinal <- function(dataname, M.var, Y.var, X.var, covariates){
  
  #Fit the gam
  gam.data <- get(dataname)
  gam.data <- gam.data %>% drop_na(all_of(c(Y.var, X.var)))
  NonNANIndex <- which(!is.na(gam.data[ ,Y.var]))
  gam.data <- gam.data[NonNANIndex,]
  tmp<-gam.data[,M.var]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp), tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  
  
  #regress covariates
  Xcorformula <- as.formula(sprintf("%s ~ %s", X.var, covariates))
  Mcorformula <- as.formula(sprintf("%s ~ %s", M.var, covariates))
  Ycorformula <- as.formula(sprintf("%s ~ %s", Y.var, covariates))
  
  gamm.model.X <- gamm4(Xcorformula, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.model.Y <- gamm4(Ycorformula, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.model.M <- gamm4(Mcorformula, random=~(1|subID), REML=TRUE, data = gam.data)
  
  
  X.res<-residuals(gamm.model.X$gam)
  M.res<-residuals(gamm.model.M$gam)
  Y.res<-residuals(gamm.model.Y$gam)
  
  X <- as.data.frame(X.res)
  Y <- as.data.frame(Y.res) 
  M <- as.data.frame(M.res)
  
  Data_Y_New <- data.frame(X=X, Y=Y, M=M)
  names(Data_Y_New) <- c("X", "Y", "M")
  
  # Set model
  model <- ' # direct effect
  Y ~ c*X
  # mediator
  M ~ a*X
  Y ~ b*M
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
  '
  
  # Run on model
  fit_sem <- sem(model, data = Data_Y_New, test="bootstrap", bootstrap=1000)
  summary(fit_sem, fit.measures=TRUE, standardize=TRUE, rsquare=TRUE)
  
  ## Calculate bootstrapped confidence intervals for the indirect (c') effect
  boot.fit <- parameterEstimates(fit_sem, boot.ci.type="perc",level=0.95, ci=TRUE,standardized = TRUE)
  #c
  direct.Beta <- boot.fit$est[1]
  direct.lwth95CI <- boot.fit$ci.lower[1]
  direct.upth95CI <- boot.fit$ci.upper[1]
  direct.P <- boot.fit$pvalue[1]
  #ab
  indirect.Beta <- boot.fit$est[7]
  indirect.Z <- boot.fit$z[7]
  indirect.lwth95CI <- boot.fit$ci.lower[7]
  indirect.upth95CI <- boot.fit$ci.upper[7]
  indirect.P <- boot.fit$pvalue[7]
  #a
  XM.Beta <- boot.fit$est[2]
  XM.lwth95CI <- boot.fit$ci.lower[2]
  XM.upth95CI <- boot.fit$ci.upper[2]
  XM.P <- boot.fit$pvalue[2]
  #b
  MY.Beta <- boot.fit$est[3]
  MY.lwth95CI <- boot.fit$ci.lower[3]
  MY.upth95CI <- boot.fit$ci.upper[3]
  MY.P <- boot.fit$pvalue[3]
  #total
  total.Beta <- boot.fit$est[8]
  total.lwth95CI <- boot.fit$ci.lower[8]
  total.upth95CI <- boot.fit$ci.upper[8]
  total.P <- boot.fit$pvalue[8]
  
  stats.results <- cbind(M.var, Y.var, direct.Beta, direct.lwth95CI, direct.upth95CI, direct.P,
                         indirect.Beta,indirect.Z, indirect.lwth95CI, indirect.upth95CI, indirect.P, 
                         XM.Beta, XM.lwth95CI, XM.upth95CI, XM.P, MY.Beta, MY.lwth95CI, MY.upth95CI, MY.P,
                         total.Beta, total.lwth95CI, total.upth95CI, total.P)
  
  return(stats.results)
}

