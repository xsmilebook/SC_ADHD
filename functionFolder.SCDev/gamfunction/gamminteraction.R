library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(lme4)
library(gamm4)
library(pbkrtest)
set.seed(925)
pboot <- function(modelobj, int_var=NA){
  numsims <- 1000
  
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- attr(terms(f1),"term.labels")
  if (!is.na(int_var)){theseVars <-  c(theseVars, int_var)}
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
  int_chisq <- pb$test["PBtest","stat"]
  
  return(c(int_pval, int_chisq))
}

#### PREDICT GAM SMOOTH FITTED VALUES FOR A SPECIFIED VALUE OF AN INTERACTING COVARIATE ####
## continuous interaction covariate
##Function to predict fitted values of a region for a given value of a covariate, using a varying coefficients smooth-by-linear covariate interaction
gamm.smooth.predict.covariateinteraction <- function(region, dataname, smooth_var, int_var, int_var.predict.percentile, covariates, knots, set_fx = FALSE, increments, stats_only=TRUE, if_residual=FALSE){
  #clean data
  gam.data <- get(dataname)
  minsmooth <- min(gam.data[,smooth_var])
  maxsmooth <- max(gam.data[,smooth_var])
  int<-gam.data[ ,int_var]
  outlierindx1<-which(int<mean(int,na.rm = T)-3*sd(int,na.rm = T) | 
                        int>mean(int,na.rm = T)+3*sd(int,na.rm = T))
  gam.data[outlierindx1 ,int_var]<-NA
  NonNANIndex <- which(!is.na(gam.data[ ,int_var]) & !is.na(gam.data[ ,region]))
  gam.data <- gam.data[NonNANIndex,]
  int_var.predict<-quantile(gam.data[ ,int_var], c(int_var.predict.percentile))
  
  # tmp<-gam.data[,region]
  # outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  # if (length(outlierindx)>0){
  #   gam.data<-gam.data[-outlierindx, ]
  # }
  parcel <- region
  
  #Fit the gam
  if (is.na(covariates)){
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s)+ s(%2$s, k=%3$s, fx=%4$s)", region, smooth_var, knots, set_fx, int_var))
    modelformula.null <- as.formula(sprintf("%1$s ~ %5$s+s(%2$s, k=%3$s, fx=%4$s)", region, smooth_var, knots, set_fx, int_var))
    modelformula.null_disease <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)", region, smooth_var, knots, set_fx))
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s)+ s(%2$s, k=%3$s, fx=%4$s) + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null <- as.formula(sprintf("%1$s ~ %5$s+s(%2$s, k=%3$s, fx=%4$s)+ %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null_disease <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)+ %5$s", region, smooth_var, knots, set_fx, covariates))
  }
  gamm.model <- gamm4(modelformula, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.results <- summary(gamm.model$gam)
  
  gamm.model.null <- gamm4(modelformula.null, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.null.results <- summary(gamm.model.null$gam)
  
  gamm.model.null_disease <- gamm4(modelformula.null_disease, random=~(1|subID), REML=TRUE, data = gam.data)
  
  T.disease <- gamm.null.results$p.table[2,3]
  P.disease <- gamm.null.results$p.table[2,4]
  
  if (stats_only==T){
    ##Full versus reduced model bootstrap p-value
    bootstrap.result.disease <- pboot(gamm.model.null)
    bootstrap.P.disease <- bootstrap.result.disease[1]
    bootstrap.chisq.disease <- bootstrap.result.disease[2]
    
    bootstrap.result.int<-pboot(gamm.model, int_var)
    bootstrap_pvalue <- bootstrap.result.int[1]
    bootstrap_chisq <- bootstrap.result.int[2]
  }else{
    bootstrap.P.disease <- NA
    bootstrap.chisq.disease <- NA
    bootstrap_pvalue<-NA
    bootstrap_chisq <- NA
  }
  gam.int.pvalue <- gamm.results$s.table[grep(x=rownames(gamm.results$s.table),pattern = int_var),"p-value"][1]
  # interaction effect size
  sse.model <- sum((gamm.model$gam$y - gamm.model$gam$fitted.values)^2)
  sse.nullmodel <- sum((gamm.model.null$gam$y - gamm.model.null$gam$fitted.values)^2)
  IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  Int.F <- gamm.results$s.table[1,3]
  
  sse.model <- sum((gamm.model.null$gam$y - gamm.model.null$gam$fitted.values)^2)
  sse.nullmodel <- sum((gamm.model.null_disease$gam$y - gamm.model.null_disease$gam$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  if(T.disease < 0){ #if the gam t-value for covariate of interest is less than 0, make the partialRsq negative
    partialRsq <- partialRsq*-1}
  
  # residuals of int_var & region
  if (is.na(covariates)){
    modelformula.int <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)", int_var, smooth_var, knots, set_fx))
    model1 <- gamm4(modelformula.int, random=~(1|subID), REML=TRUE, data = gam.data)
    res.int_var <- residuals(model1$gam)
    modelformula.region <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)", region, smooth_var, knots, set_fx))
    model2 <- gamm4(modelformula.region, random=~(1|subID), REML=TRUE, data = gam.data)
    res.region <- residuals(model2$gam)
  }else{
    modelformula.int <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)+ %5$s", int_var, smooth_var, knots, set_fx, covariates))
    model1 <- gamm4(modelformula.int, random=~(1|subID), REML=TRUE, data = gam.data)
    res.int_var <- residuals(model1$gam)
    modelformula.region <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)+ %5$s", region, smooth_var, knots, set_fx, covariates))
    model2 <- gamm4(modelformula.region, random=~(1|subID), REML=TRUE, data = gam.data)
    res.region <- residuals(model2$gam)
  }

  ### get predicted data
  ##############################
  modelobj <- gamm.model$gam
  #Extract gam input data
  df <- modelobj$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  np <- increments #number of predictions to make; predict at np increments of smooth_var
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(modelobj$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(modelobj$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(modelobj$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(minsmooth,maxsmooth, length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on the modal number
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on the modal number
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init)
  pred[,int_var] <- as.numeric(int_var.predict)
  
  #Generate fitted (predicted) values based on the gam model and predication data frame
  predicted.smooth <- fitted_values(object = modelobj, data = pred)
  predicted.smooth$fitted.centered <- scale(predicted.smooth$.fitted, center=T, scale = F) #subtract the intercept from fitted values
  predicted.smooth <- predicted.smooth %>% dplyr::select(all_of(smooth_var), .fitted, .se, .lower_ci, .upper_ci, fitted.centered)
  
  changed.range <- predicted.smooth$.fitted[which.max(predicted.smooth$age)]-predicted.smooth$.fitted[which.min(predicted.smooth$age)]
  changed.ratio <- predicted.smooth$.fitted[which.max(predicted.smooth$age)] / predicted.smooth$.fitted[which.min(predicted.smooth$age)]
  SCweight.agemax <- predicted.smooth$.fitted[which.max(predicted.smooth$age)]
  SCweight.agemin <- predicted.smooth$.fitted[which.min(predicted.smooth$age)]
  ############################################
  
  stats.reults<-cbind(parcel, int_var, bootstrap_pvalue, bootstrap_chisq, gam.int.pvalue, IntpartialRsq,partialRsq, Int.F, 
                      T.disease, P.disease, bootstrap.P.disease, bootstrap.chisq.disease, changed.range, 
                      changed.ratio, SCweight.agemax, SCweight.agemin)
  
  full.results<-list(stats.reults, predicted.smooth)
  if(if_residual == TRUE)
    return(data.frame(region.res=res.region, int_var.res=res.int_var, scanID=gam.data$scanID))
  if(stats_only == TRUE)
    return(stats.reults)
  if(stats_only == FALSE)
    return(full.results)
  
}