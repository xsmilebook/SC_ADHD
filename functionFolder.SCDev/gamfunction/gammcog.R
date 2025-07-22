library(tidyr)
library(mgcv)
library(psych)
library(tidyverse)
library(dplyr)
library(lme4)
library(gamm4)
library(pbkrtest)

pboot <- function(modelobj){
  numsims <- 1000
  
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- attr(terms(f1),"term.labels")
  f2 <- reformulate(theseVars[2:(length(theseVars))],response = thisResp)
  
  g1 <- gam(f1,data = df)
  g2 <- gam(f2,data = df)
  
  mat1 <- model.matrix(g1)
  mat2 <- model.matrix(g2)
  
  subID<-df$subID
  y <- df[,thisResp]
  
  m1 <- lmer(y ~ -1 + mat1 + (1|subID))
  m2 <- lmer(y ~ -1 + mat2 + (1|subID))
  refdist <- PBrefdist(m1, m2, nsim=numsims)#, cl=cl)
  pb <- PBmodcomp(m1, m2, ref = refdist)
  int_pval <- pb$test["PBtest","p.value"]
  
  return(int_pval)
}



#### FIT GAM SMOOTH ####
##Function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
gam.fit.cognition <- function(region, dataname, cognition_var, smooth_var, covariates, knots,corrmethod, set_fx = FALSE, stats_only = FALSE){
  
  #Fit the gam
  gam.data <- get(dataname)
  cognition<-gam.data[ ,cognition_var]
  outlierindx1<-which(cognition<mean(cognition,na.rm = T)-3*sd(cognition,na.rm = T) | cognition>mean(cognition,na.rm = T)+3*sd(cognition,na.rm = T))
  gam.data[outlierindx1 ,cognition_var]<-NA
  NonNANIndex <- which(!is.na(gam.data[ ,cognition_var]) & !is.na(gam.data[ ,region]))
  
  gam.data <- gam.data[NonNANIndex,]
  tmp<-gam.data[,region]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  parcel <- as.character(region)
  modelformula <- as.formula(sprintf("%s ~ %s + s(%s, k = %s, fx = %s) + %s", region, cognition_var, smooth_var, knots, set_fx, covariates))
  modelformula.null<-as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s",region, smooth_var, knots, set_fx, covariates))
  gamm.model <- gamm4(modelformula, random=~(1|subID), REML=FALSE, data = gam.data)
  gamm.results <- summary(gamm.model$gam)
  gamm.model.null <- gamm4(modelformula.null, random=~(1|subID), REML=FALSE, data = gam.data)
  
  #Send to bootstrap function
  bootstrap_pvalue<-pboot(gamm.model) 

  #GAM statistics
  #F value for the smooth term and GAM-based significance of the smooth term
  gamm.smooth.t <- gamm.results$p.table[2,3]
  gamm.smooth.pvalue <- gamm.results$p.table[2,4]
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gamm.model$gam$y - gamm.model$gam$fitted.values)^2)
  sse.nullmodel <- sum((gamm.model.null$gam$y - gamm.model.null$gam$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  ### effect direction
  if(gamm.smooth.t < 0){ #if the gam t-value for covariate of interest is less than 0, make the partialRsq negative
    partialRsq <- partialRsq*-1
    }
  
  #residual correlation
  varcorformula1 <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  varcorformula2 <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", cognition_var, smooth_var, knots, set_fx, covariates))
  res1<-residuals(gamm4(varcorformula1, random=~(1|subID), REML=FALSE, data = gam.data)$gam)
  res2<-residuals(gamm4(varcorformula2, random=~(1|subID), REML=FALSE, data = gam.data)$gam)
  PCorr_Test <- corr.test(res1, res2, method=corrmethod)
  correstimate<-PCorr_Test$r
  corrp <- PCorr_Test$p
  
  stats.results <- cbind(parcel, cognition_var, gamm.smooth.t, gamm.smooth.pvalue,bootstrap_pvalue, partialRsq, correstimate, corrp)
  data.results <- list()
  data.results[[1]] <- as.data.frame(stats.results)
  data.results[[2]] <- data.frame(SCres=res1, cogres=res2)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(data.results)
}