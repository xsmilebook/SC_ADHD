library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)

#### PREDICT GAM SMOOTH FITTED VALUES FOR A SPECIFIED VALUE OF AN INTERACTING COVARIATE ####
##Function to predict fitted values of a region for a given value of a covariate, using a varying coefficients smooth-by-linear covariate interaction
gam.factor.predict.covariateinteraction <- function(region, dataname, smooth_var, int_var,cog_var, covariates, knots, set_fx = FALSE, increments, mod_need=FALSE){
  #Fit the gam
  gam.data <- get(dataname)
  tmp<-gam.data[,region]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp), tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  
  gam.data <- get(dataname)
  tmp<-gam.data[,cog_var]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp), tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  
  parcel <- region
  
  #Fit the gam
  modelformula <- as.formula(sprintf("%1$s ~ %5$s*%6$s +s(%2$s, k=%3$s, fx=%4$s)+ %7$s", region, smooth_var, knots, set_fx, int_var, cog_var, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)+%5$s+ %6$s+ %7$s", region, smooth_var, knots, set_fx, int_var, cog_var, covariates))
  gam.model.null <- gam(modelformula.null, method = "REML", data = gam.data)
  gam.null.results <- summary(gam.model.null)
  
  ##Full versus reduced model anova p-value
  anova.int.pvalue <- anova(gam.model.null, gam.model,test='Chisq')$`Pr(>Chi)`[2]
  
  #GAM derivatives
  levels<-levels(gam.data$gender)
  gam.data.F <-gam.data[gam.data$gender==levels[1], ]
  gam.data.M <-gam.data[gam.data$gender==levels[2], ]
  modelformula.sep <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)+ %5$s", region, smooth_var, knots, set_fx, covariates))
  gam.model.F <- gam(modelformula.sep, method = "REML", data = gam.data.F)
  gam.model.M <- gam(modelformula.sep, method = "REML", data = gam.data.M)
  
  #Get derivatives of the smooth function using finite differences
  derv.F <- derivatives(gam.model.F, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = F) #derivative at 200 indices of smooth_var with a simultaneous CI
  #Identify derivative significance window(s)
  derv.F <- derv.F %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv.F$sig_deriv = derv.F$derivative*derv.F$sig
  derv.M <- derivatives(gam.model.M, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = F) #derivative at 200 indices of smooth_var with a simultaneous CI
  derv.M <- derv.M %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv.M$sig_deriv = derv.M$derivative*derv.M$sig
  
  if (!is.na(anova.int.pvalue) & anova.int.pvalue<0.05){
    ## effect size
    sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
    sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
    IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
    anova.int.zvalue = qnorm(anova.int.pvalue / 2, lower.tail=FALSE)
    edf.F <-gam.results$s.table[1,1]
    edf.M <-gam.results$s.table[2,1]
    gam.pvalue.F<-gam.results$s.table[1,4]
    gam.pvalue.M<-gam.results$s.table[2,4]
    gam.Fvalue.F<-gam.results$s.table[1,3]
    gam.Fvalue.M<-gam.results$s.table[2,3]
    #### F
    #Age of decrease onset and offset
    if(sum(derv.F$sig) > 0){ 
      decreasing.range.F <- derv.F$data[derv.F$sig_deriv < 0] #identify all ages with a significant negative derivative (i.e., smooth_var indices where y is decreasing)
      if(length(decreasing.range.F) > 0){
        decrease.onset.F <- min(decreasing.range.F) #find youngest age with a significant negative derivative
        decrease.offset.F <- max(decreasing.range.F)
      }
      if(length(decreasing.range.F) == 0){
        decrease.onset.F <- NA
        decrease.offset.F <-NA
      }}
    if(sum(derv.F$sig) == 0){
      decrease.onset.F <- NA
      decrease.offset.F <-NA}  
    
    #Age of increase onset and offset
    if(sum(derv.F$sig) > 0){ 
      increasing.range.F <- derv.F$data[derv.F$sig_deriv > 0] #identify all ages with a significant positive derivative (i.e., smooth_var indices where y is increasing)
      if(length(increasing.range.F) > 0){
        increase.onset.F <- min(increasing.range.F)
        increase.offset.F <- max(increasing.range.F) #find oldest age with a significant positive derivative
        
      }
      if(length(increasing.range.F) == 0){
        increase.onset.F <- NA
        increase.offset.F <- NA
      }
    }
    if(sum(derv.F$sig) == 0){ 
      increase.onset.F <- NA
      increase.offset.F <- NA}
    ### M
    #Age of decrease onset and offset
    if(sum(derv.M$sig) > 0){ 
      decreasing.range.M <- derv.M$data[derv.M$sig_deriv < 0] #identify all ages with a significant negative derivative (i.e., smooth_var indices where y is decreasing)
      if(length(decreasing.range.M) > 0){
        decrease.onset.M <- min(decreasing.range.M) #find youngest age with a significant negative derivative
        decrease.offset.M <- max(decreasing.range.M)
      }
      if(length(decreasing.range.M) == 0){
        decrease.onset.M <- NA
        decrease.offset.M <-NA
      }}
    if(sum(derv.M$sig) == 0){
      decrease.onset.M <- NA
      decrease.offset.M <-NA}  
    
    #Age of increase onset and offset
    if(sum(derv.M$sig) > 0){ 
      increasing.range.M <- derv.M$data[derv.M$sig_deriv > 0] #identify all ages with a significant positive derivative (i.e., smooth_var indices where y is increasing)
      if(length(increasing.range.M) > 0){
        increase.onset.M <- min(increasing.range.M)
        increase.offset.M <- max(increasing.range.M) #find oldest age with a significant positive derivative
        
      }
      if(length(increasing.range.M) == 0){
        increase.onset.M <- NA
        increase.offset.M <- NA
      }
    }
    if(sum(derv.M$sig) == 0){ 
      increase.onset.M <- NA
      increase.offset.M <- NA}
    
  }else{
    IntpartialRsq<-NA
    anova.int.zvalue = NA
    edf.F <-NA
    edf.M <-NA
    gam.pvalue.F<-NA
    gam.pvalue.M<-NA
    gam.Fvalue.F<-NA
    gam.Fvalue.M<-NA
    decrease.onset.F <- NA
    decrease.offset.F <-NA
    decrease.onset.M <- NA
    decrease.offset.M <-NA
    increase.onset.F <- NA
    increase.offset.F <- NA
    increase.onset.M <- NA
    increase.offset.M <- NA
  }
  
  full.results<-cbind(parcel, anova.int.pvalue, IntpartialRsq, anova.int.zvalue, edf.F, edf.M, gam.pvalue.F, gam.pvalue.M,
                      gam.Fvalue.F, gam.Fvalue.M, decrease.onset.F, decrease.onset.M, increase.onset.F, increase.onset.M,
                      decrease.offset.F, decrease.offset.M, increase.offset.F, increase.offset.M)
  gam.model<-list()
  gam.model[[1]]<-gam.data.F
  gam.model[[2]]<-gam.data.M
  if(mod_need == TRUE)
    return(gam.model)
  if(mod_need == FALSE)
    return(full.results)
}