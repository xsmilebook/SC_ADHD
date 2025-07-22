library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)


#### FIT GAM SMOOTH ####
##Function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx))) per each region in atlas and save out statistics and derivative-based characteristics
gam.fit.smooth <- function(depvar, region, dataname, smooth_var, knots, set_fx = FALSE, stats_only = FALSE, mod_only = FALSE){
  
  #Fit the gam
  gam.data <- get(dataname)
  gam.data <- gam.data[[region]]
  #parcel <- paste0("SC.", as.character(region))
  parcel <- as.character(region)
  #depvar <- sprintf("%s.%s", 'SC', as.character(region))
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s)", depvar, smooth_var, knots, set_fx))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #GAM derivatives
  #Get derivatives of the smooth function using finite differences
  derv <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = F) #derivative at 200 indices of smooth_var with a simultaneous CI
  derv2 <- derivatives(gam.model,order = 2L, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = F)
  derv2 <- derv2 %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) 
  derv2$sig_deriv = derv2$derivative*derv2$sig
  meanderv2<-mean(derv2$sig_deriv)
  #Identify derivative significance window(s
  derv <- derv %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv$sig_deriv = derv$derivative*derv$sig #add "sig_deriv derivatives column where non-significant derivatives are set to 0
  
  #GAM statistics
  #F value for the smooth term and GAM-based significance of the smooth term
  gam.smooth.F <- gam.results$s.table[3]
  gam.smooth.pvalue <- gam.results$s.table[4]
  gam.edf <- gam.results$s.table[1]
  
  ##Full versus reduced model direction-dependent partial R squared
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  Rsq <- gam.results$r.sq
  
  res1<-gam.data[,depvar]
  res2<-gam.data[,smooth_var]
  PCorr_Test <- cor.test(res1, res2)
  correstimate<-PCorr_Test$estimate
  corrp <- PCorr_Test$p.value
  ### effect direction
  mean.derivative <- mean(derv$derivative)
  if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
    Rsq <- Rsq*-1}
  
  #Derivative-based temporal characteristics
  #Age of developmental change onset
  if(sum(derv$sig) > 0){ #if derivative is significant at at least 1 age
    change.onset <- min(derv$data[derv$sig==T])} #find first age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ #if gam derivative is never significant
    change.onset <- NA} #assign NA
  
  #Age of maximal developmental change
  if(sum(derv$sig) > 0){ 
    derv$abs_sig_deriv = round(abs(derv$sig_deriv),5) #absolute value significant derivatives
    maxval <- max(derv$abs_sig_deriv) #find the largest derivative
    window.peak.change <- derv$data[derv$abs_sig_deriv == maxval] #identify the age(s) at which the derivative is greatest in absolute magnitude
    peak.change <- mean(window.peak.change)} #identify the age of peak developmental change
  if(sum(derv$sig) == 0){ 
    peak.change <- NA}  
  
  #Age of maximal increasing change
  if(sum(derv$sig) > 0){
    derv$sig_deriv_round<-round(derv$sig_deriv, 5)
    increasing.derv <- derv$sig_deriv_round[derv$sig_deriv > 0]
    if(length(increasing.derv) > 0){
      maxval.increase <- max(increasing.derv)
      window.peak.increase.change <- derv$data[derv$sig_deriv_round == maxval.increase]
      peak.increase.change <- mean(window.peak.increase.change)
    }else{peak.increase.change<-NA}}else{peak.increase.change<-NA}
  
  #Age of decrease onset and offset
  if(sum(derv$sig) > 0){ 
    decreasing.range <- derv$data[derv$sig_deriv < 0] #identify all ages with a significant negative derivative (i.e., smooth_var indices where y is decreasing)
    if(length(decreasing.range) > 0){
      decrease.onset <- min(decreasing.range) #find youngest age with a significant negative derivative
      decrease.offset <- max(decreasing.range)}
    if(length(decreasing.range) == 0){
      decrease.onset <- NA
      decrease.offset <-NA
    }}
  if(sum(derv$sig) == 0){
    decrease.onset <- NA
    decrease.offset <-NA}  
  
  
  #Age of increase onset and offset
  if(sum(derv$sig) > 0){ 
    increasing.range <- derv$data[derv$sig_deriv > 0] #identify all ages with a significant positive derivative (i.e., smooth_var indices where y is increasing)
    if(length(increasing.range) > 0){
      increase.onset <- min(increasing.range)
      increase.offset <- max(increasing.range) #find oldest age with a significant positive derivative
      
    }
    if(length(increasing.range) == 0){
      increase.onset <- NA
      increase.offset <- NA
    }
  }
  if(sum(derv$sig) == 0){ 
    increase.onset <- NA
    increase.offset <- NA}  
  
  #Age of maturation
  if(sum(derv$sig) > 0){ 
    change.offset <- max(derv$data[derv$sig==T])} #find last age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ 
    change.offset <- NA}  
  
  full.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, Rsq,correstimate, 
                        corrp,gam.edf, meanderv2, change.onset, peak.change,peak.increase.change, 
                        decrease.onset, decrease.offset, increase.onset, 
                        increase.offset, change.offset)
  stats.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, Rsq, correstimate, corrp,gam.edf, meanderv2)
  if(mod_only == TRUE)
    return(gam.model)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(full.results)
}