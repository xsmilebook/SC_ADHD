library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)
#### CALCULATE SMOOTH ESTIMATES ####
##Function to estimate the zero-averaged gam smooth function 
gam.estimate.smooth <- function(dataname, region, smooth_var, covariates, knots, set_fx = FALSE, increments){
  
  #Fit the gam
  gam.data <- get(dataname)
  parcel <- region
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  np <- increments #number of predictions to make; predict at np increments of smooth_var
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]},
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]}
            )
    }
  }
  pred <- thisPred %>% dplyr::select(-init)
  
  #Estimate the smooth trajectory 
  estimated.smooth <- smooth_estimates(object = gam.model, data = pred)
  estimated.smooth <- estimated.smooth %>% dplyr::select(age, est)
  
  return(estimated.smooth)
}

