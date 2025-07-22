gamm.varyingcoefficients <- function(region, dataname, smooth_var, int_var, covariates, knots, set_fx = FALSE, increments, draws, return_posterior_coefficients = FALSE){
  
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  #Fit the gam
  gam.data <- get(dataname)
  smoothmin <- min(gam.data[ , smooth_var])
  smoothmax <- max(gam.data[ , smooth_var])
  cognition<-gam.data[ ,int_var]
  outlierindx1<-which(cognition<mean(cognition,na.rm = T)-3*sd(cognition,na.rm = T) | cognition>mean(cognition,na.rm = T)+3*sd(cognition,na.rm = T))
  gam.data[outlierindx1 ,int_var]<-NA
  NonNANIndex <- which(!is.na(gam.data[ ,int_var]) & !is.na(gam.data[ ,region]))
  
  gam.data <- gam.data[NonNANIndex,]
  tmp<-gam.data[,region]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  parcel <- region
  modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + s(%2$s, by=%5$s, k=%3$s, fx=%4$s) + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  gamm.model <- gamm4(modelformula, random=~(1|subID), REML=TRUE, data = gam.data)
  
  #Extract gam input data
  df <- gamm.model$gam$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame, used to estimate (posterior) model slopes (varying covariate coefficients)
  theseVars <- attr(gamm.model$gam$terms,"term.labels") 
  varClasses <- attr(gamm.model$gam$terms,"dataClasses") 
  
  #prediction df for int_var P10th at np smooth_var increments
  pred.low <- data.frame(init = rep(0,np)) 
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == int_var) { 
      #pred.low[,int_var] <- (min(df[,int_var],na.rm = T))
      pred.low[,int_var] <- (quantile(df[,int_var], probs=0.1, na.rm = T))
    } else if (thisVar == smooth_var) {
      pred.low[,smooth_var] = seq(smoothmin,smoothmax, length.out = np)
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {pred.low[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {pred.low[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on first level of factor 
              "ordered" = {pred.low[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on first level of ordinal variable
      )}}
  
  #prediction df for int_var P90th at np smooth_var increments
  pred.high <- data.frame(init = rep(0,np)) 
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == int_var) { 
      #pred.high[,int_var] <- (max(df[,int_var],na.rm = T))
      pred.high[,int_var] <- (quantile(df[,int_var], probs=0.9, na.rm = T))
    } else if (thisVar == smooth_var) {
      pred.high[,smooth_var] = seq(smoothmin,smoothmax, length.out = np)
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {pred.high[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {pred.high[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on first level of factor 
              "ordered" = {pred.high[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on first level of ordinal variable
      )}}
  
  pred <- rbind(pred.low, pred.high) #complete pred df 
  pred <- pred %>% dplyr::select(-init)
  
  #Get effects (slopes) along the smooth function for the true model
  if(return_posterior_coefficients == FALSE){
    #varying coefficient slopes
    predicted.values <- fitted_values(object = gamm.model$gam, data = pred) #predict y at min and mix int_var along the smooth function
    predicted.values <- predicted.values %>% dplyr::select(fitted, all_of(smooth_var), all_of(int_var))
    colnames(predicted.values) <- c("fitted", "smooth.var", "int.var")
    predicted.values$smooth.var <- round(predicted.values$smooth.var, 3)
    
    varyingcoeff.slopes <-  predicted.values %>% #calculate the effect of int_var on y  (slope; delta y/delta int_var) along the smooth function 
      group_by(smooth.var) %>%
      do(slope = diff(.$fitted)/diff(.$int.var)) %>%
      unnest(cols = c(slope))
    colnames(varyingcoeff.slopes) <- c(smooth_var, "slope")  
  }
  
  #Estimate posterior distribution of effects (slopes) from simulated GAM posterior distribution
  if(return_posterior_coefficients == TRUE){
    Vb <- vcov(gamm.model$gam, unconditional = UNCONDITIONAL) #variance-covariance matrix of fitted gam coefficients
    sims <- MASS::mvrnorm(npd, mu = coef(gamm.model$gam), Sigma = Vb) #simulated model coefficients for npd draws from the posterior
    X0 <- predict(gamm.model$gam, newdata = pred, type = "lpmatrix") #get matrix of linear predictors that maps model parameters to the smooth fit (outcome measure scale)
    predicted.values.posterior <- X0 %*% t(sims) #predicted/fitted values along smooth_var, i.e., posterior smooths
    
    predicted.values.posterior <- as.data.frame(predicted.values.posterior)
    colnames(predicted.values.posterior) <- sprintf("draw%s",seq(from = 1, to = npd))
    predicted.values.posterior <- cbind(as.numeric(pred[,smooth_var]), as.numeric(pred[,int_var]), predicted.values.posterior)
    colnames(predicted.values.posterior)[1] <- c("smooth.var")
    colnames(predicted.values.posterior)[2] <- c("int.var")
    predicted.values.posterior$smooth.var <- round(predicted.values.posterior$smooth.var, 3)
    
    varyingcoeff.slopes.CI = predicted.values.posterior %>% pivot_longer(cols = contains("draw"), names_to = "draw",values_to = "posterior.fitted") 
    varyingcoeff.slopes.CI$int.var[varyingcoeff.slopes.CI$int.var == quantile(df[,int_var], 0.1)] <- c("low")
    varyingcoeff.slopes.CI$int.var[varyingcoeff.slopes.CI$int.var == quantile(df[,int_var], 0.9)] <- c("high")
    varyingcoeff.slopes.CI <- varyingcoeff.slopes.CI %>% pivot_wider(names_from = "int.var", values_from = "posterior.fitted") %>% mutate(slope = (high-low)/(quantile(df[,int_var], 0.9) - quantile(df[,int_var], 0.1))) 
    #calculate the effect of int_var on y  (slope; delta y/delta int_var) along the smooth function for all draws
    
    varyingcoeff.slopes.CI <- varyingcoeff.slopes.CI %>% dplyr::select(smooth.var, draw, slope)
    varyingcoeff.slopes.CI <- cbind(as.character(parcel), varyingcoeff.slopes.CI)
    colnames(varyingcoeff.slopes.CI) <- c("label", smooth_var, "draw", "slope")
  }
  
  if(return_posterior_coefficients == FALSE)
    return(varyingcoeff.slopes)
  if(return_posterior_coefficients == TRUE)
    return(varyingcoeff.slopes.CI)

}