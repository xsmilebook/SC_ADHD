library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(ecostats)

#### Interaction analysis for linear effects (cognition) by  a categorical variable####
##Function to predict fitted values of a region for a each level of a categorical variable, using a varying coefficients linear-by-category interaction
gam.smooth.predict.interaction <- function(region, dataname, smooth_var, int_var, covariates, knots, set_fx = FALSE, stats_only=FALSE, increments=1000){
  #Fit the gam
  gam.data <- get(dataname)
  
  parcel <- region
  
  if (is.na(covariates)){
    modelformula.main.effect.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)", region, smooth_var, knots, set_fx))
  } else {
    modelformula.main.effect.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + %5$s", region, smooth_var, knots, set_fx, covariates))
  }
  
  #Fit the gam
  if (is.na(covariates)){
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s)+ %5$s", region, smooth_var, knots, set_fx, int_var))
    modelformula.null <- as.formula(sprintf("%1$s ~ %5$s+s(%2$s, k=%3$s, fx=%4$s)", region, smooth_var, knots, set_fx, int_var))
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s) + %5$s + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null <- as.formula(sprintf("%1$s ~ %5$s+s(%2$s, k=%3$s, fx=%4$s)+ %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  }

  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  gam.model.null <- gam(modelformula.null, method = "REML", data = gam.data)
  gam.null.results <- summary(gam.model.null)
  
  gam.model.main.effect.null <- gam(modelformula.main.effect.null, method="REML", data=gam.data)
  
  ##Full versus reduced model anova p-value
  #anova.int.pvalue <- anova(gam.model.null, gam.model,test='Chisq')$`Pr(>Chi)`[2]
  cat("\n--- Running anovaPB for parcel:", region, "---\n")
  intresults <- try(anovaPB(gam.model.null, gam.model, n.sim = 1000,test='Chisq', ncpus=1))
  cat("anovaPB results:\n")
  print(intresults)
  
  anova.int.pvalue <- intresults$"Pr(>Chi)"[2]
  anova.int.chisq <- intresults$"Deviance"[2]
  
  # interaction effect pvalue
  gam.int.pvalue <- gam.results$s.table[grep(x=rownames(gam.results$s.table),pattern = int_var),"p-value"]
  # interaction effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  
  p_table_null <- gam.null.results$p.table
  
  factor_levels <- levels(gam.data[[int_var]])
  non_reference_level <- factor_levels[2] 
  
  main_effect_row_name <- paste0(int_var, non_reference_level)
  T.disease <- p_table_null[main_effect_row_name, "t value"]
  P.disease <- p_table_null[main_effect_row_name, "Pr(>|t|)"]

  
  cat("\n--- Running anovaPB for MAIN EFFECT in parcel:", region, "---\n")
  maineffect_results <- try(anovaPB(gam.model.main.effect.null, gam.model.null, n.sim=1000, ncpus=1))
  print(maineffect_results)
  if (inherits(maineffect_results, "try-error")) {
    bootstrap.P.disease <- NA
    bootstrap.chisq.disease <- NA
  } else {
    bootstrap.P.disease <- maineffect_results$"Pr(>F)"[2]
    bootstrap.chisq.disease <- maineffect_results$"Deviance"[2]
  }
  
  ### get predicted data
  ##############################
  modelobj <- gam.model
  #Extract gam input data
  df <- modelobj$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  np <- increments #number of predictions to make; predict at np increments of smooth_var
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(modelobj$terms,"term.labels") #gam model predictors (smooth_var + covariates)

  varClasses <- attr(modelobj$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(modelobj$terms[[2]]) #the measure to predict
  for (thisVar in theseVars) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "nmatrix.1" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on the modal number
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on the modal number
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init)
  
  original_levels <- levels(modelobj$model[[int_var]])
  level_ref <- original_levels[1]
  level_comp <- original_levels[2]
  pred.HC <- pred.Disease <- pred
  pred.HC[[int_var]] <- factor(level_ref, levels = original_levels)
  pred.Disease[[int_var]] <- factor(level_comp, levels = original_levels)
  
  #Generate fitted (predicted) values based on the gam model and predication data frame
  predicted.smooth.HC <- fitted_values(object = modelobj, data = pred.HC)

  if (is.null(predicted.smooth.HC) || !".fitted" %in% names(predicted.smooth.HC) || is.null(predicted.smooth.HC$.fitted)) {
    stop("fitted_values() failed to return a valid 'fitted' column for HC.")
  }
  predicted.smooth.HC$fitted.centered <- scale(predicted.smooth.HC$.fitted, center=T, scale = F) #subtract the intercept from fitted values
  predicted.smooth.Disease <- fitted_values(object = modelobj, data = pred.Disease)
  predicted.smooth.Disease$fitted.centered <- scale(predicted.smooth.Disease$.fitted, center=T, scale = F) #subtract the intercept from fitted values
  predicted.smooth <- rbind(predicted.smooth.HC, predicted.smooth.Disease)
  

  stats.reults <- data.frame(
    parcel = parcel,
    int_var = int_var,
    anova.int.pvalue = anova.int.pvalue,
    anova.int.chisq = anova.int.chisq,  
    IntpartialRsq = IntpartialRsq,          
    T.disease = T.disease,               
    P.disease = P.disease,           
    bootstrap.P.disease = bootstrap.P.disease, 
    bootstrap.chisq_disease = bootstrap.chisq.disease
  )
  full.results<-list(stats.reults, predicted.smooth)
  if(stats_only == TRUE)
    return(stats.reults)
  if(stats_only == FALSE)
    return(full.results)
}