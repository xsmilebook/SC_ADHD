## This function compute the effect on age derivative at each specific age point.
#### DERIVATIVES ####
library(gratia)

##Function to compute smooth derivatives for a main GAM model and for individual draws from the simulated posterior distribution
gam.derivatives.changed <- function(region, dataname, smooth_var,int_var,covariates, knots, set_fx = FALSE, draws, increments, return_posterior_derivatives = TRUE){
  
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior derivative sets estimated
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  EPS <- 1e-07 #finite differences
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  #Fit the gam
  gam.data <- get(dataname)
  parcel <- region
  modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + s(%2$s, by=%5$s, k=%3$s, fx=%4$s) + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  gam.model <- gam(modelformula, method = "REML", data=gam.data)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  theseVars <- attr(gam.model$terms,"term.labels")
  varClasses <- attr(gam.model$terms,"dataClasses")
  thisResp <- as.character(gam.model$terms[[2]])
  #Create a prediction data frame, used to estimate (posterior) model coefficients
  #### low level
  thisPred.low <- data.frame(init = rep(0,np)) 
  
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) {
      thisPred.low[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
    } else if (thisVar == int_var) {
      thisPred.low[,int_var] <- (quantile(df[,int_var], probs=0.1, na.rm = T))
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred.low[,thisVar] = median(df[,thisVar])},
              "factor" = {thisPred.low[,thisVar] = levels(df[,thisVar])[[levelact]]},
              "ordered" = {thisPred.low[,thisVar] = levels(df[,thisVar])[[levelact]]}
      )
    }
  }
  pred.low <- thisPred.low %>% dplyr::select(-init)
  pred2.low <- pred.low #second prediction df
  pred2.low[,smooth_var] <- pred.low[,smooth_var] + EPS #finite differences
  
  #### high level
  thisPred.high <- data.frame(init = rep(0,np)) 
  
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) {
      thisPred.high[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
    } else if (thisVar == int_var) {
      thisPred.high[,int_var] <- (quantile(df[,int_var], probs=0.9, na.rm = T))
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred.high[,thisVar] = median(df[,thisVar])},
              "factor" = {thisPred.high[,thisVar] = levels(df[,thisVar])[[levelact]]},
              "ordered" = {thisPred.high[,thisVar] = levels(df[,thisVar])[[levelact]]}
      )
    }
  }
  pred.high <- thisPred.high %>% dplyr::select(-init)
  pred2.high <- pred.high #second prediction df
  pred2.high[,smooth_var] <- pred.high[,smooth_var] + EPS #finite differences
  
  #Estimate smooth derivatives
  derivs.low <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = UNCONDITIONAL, data = pred.low) #derivative at 200 indices of smooth_var with a simultaneous CI
  derivs.fulldf.low <- derivs.low %>% select(data, derivative, se, lower, upper)
  derivs.fulldf.low <- derivs.fulldf.low %>% mutate(significant = !(0 > lower & 0 < upper))
  derivs.fulldf.low$significant.derivative = derivs.fulldf.low$derivative*derivs.fulldf.low$significant
  colnames(derivs.fulldf.low) <- c(sprintf("%s", smooth_var), "derivative", "se", "lower", "upper", "significant", "significant.derivative")
  
  derivs.high <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = UNCONDITIONAL, data = pred.high) #derivative at 200 indices of smooth_var with a simultaneous CI
  derivs.fulldf.high <- derivs.high %>% select(data, derivative, se, lower, upper)
  derivs.fulldf.high <- derivs.fulldf.high %>% mutate(significant = !(0 > lower & 0 < upper))
  derivs.fulldf.high$significant.derivative = derivs.fulldf.high$derivative*derivs.fulldf.high$significant
  colnames(derivs.fulldf.high) <- c(sprintf("%s", smooth_var), "derivative", "se", "lower", "upper", "significant", "significant.derivative")
  
  derivs.fulldf.diff <- derivs.fulldf.low[,c(1)]
  derivs.fulldf.diff$derivative.changed <- derivs.fulldf.high$derivative - derivs.fulldf.low$derivative
  #Estimate posterior smooth derivatives from simulated GAM posterior distribution
  if(return_posterior_derivatives == TRUE){
    Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
    sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance
    ## low level
    X0 <- predict(gam.model, newdata = pred.low, type = "lpmatrix") #get matrix of linear predictors for pred
    X1 <- predict(gam.model, newdata = pred2.low, type = "lpmatrix") #get matrix of linear predictors for pred2
    Xp <- (X1 - X0) / EPS 
    posterior.derivs.low <- Xp %*% t(sims) #Xp * simulated model coefficients = simulated derivatives. Each column of posterior.derivs contains derivatives for a different draw from the simulated posterior distribution
    posterior.derivs.low <- as.data.frame(posterior.derivs.low)
    # posterior.derivs[1,1]<-sum(Xp[1,] * sims[1,]), posterior.derivs[1,2]<-sum(Xp[1,] * sims[2,])...
    # posterior.derivs is the change of posterior smooth term when there is an infinitesimal change of age
    # namely posterior derivatives
    colnames(posterior.derivs.low) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
    posterior.derivs.low <- cbind(as.numeric(pred.low[,smooth_var]), posterior.derivs.low) #add smooth_var increments from pred df to first column
    colnames(posterior.derivs.low)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
    posterior.derivs.low <- cbind(as.character(region), posterior.derivs.low) #add parcel label to first column
    colnames(posterior.derivs.low)[1] <- "label_ID" #label the column
    posterior.derivs.low.long <- posterior.derivs.low %>% pivot_longer(contains("draw"), names_to = "draw",values_to = "posterior.derivative")
    
    ## high level
    X0 <- predict(gam.model, newdata = pred.high, type = "lpmatrix") #get matrix of linear predictors for pred
    X1 <- predict(gam.model, newdata = pred2.high, type = "lpmatrix") #get matrix of linear predictors for pred2
    Xp <- (X1 - X0) / EPS 
    posterior.derivs.high <- Xp %*% t(sims) #Xp * simulated model coefficients = simulated derivatives. Each column of posterior.derivs contains derivatives for a different draw from the simulated posterior distribution
    posterior.derivs.high <- as.data.frame(posterior.derivs.high)
    # posterior.derivs[1,1]<-sum(Xp[1,] * sims[1,]), posterior.derivs[1,2]<-sum(Xp[1,] * sims[2,])...
    # posterior.derivs is the change of posterior smooth term when there is an infinitesimal change of age
    # namely posterior derivatives
    colnames(posterior.derivs.high) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
    posterior.derivs.high <- cbind(as.numeric(pred.high[,smooth_var]), posterior.derivs.high) #add smooth_var increments from pred df to first column
    colnames(posterior.derivs.high)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
    posterior.derivs.high <- cbind(as.character(region), posterior.derivs.high) #add parcel label to first column
    colnames(posterior.derivs.high)[1] <- "label_ID" #label the column
    posterior.derivs.high.long <- posterior.derivs.high %>% pivot_longer(contains("draw"), names_to = "draw",values_to = "posterior.derivative")
    
    posterior.derivs.diff.long <- posterior.derivs.low.long[,c(1:3)]
    posterior.derivs.diff.long$posterior.derivs.changed <- 
      posterior.derivs.high.long$posterior.derivative-posterior.derivs.low.long$posterior.derivative
  } #np*npd rows, 3 columns (smooth_var, draw, posterior.derivative)
  
  if(return_posterior_derivatives == FALSE)
    return(derivs.fulldf.diff)
  if(return_posterior_derivatives == TRUE)
    return(posterior.derivs.diff.long)
}