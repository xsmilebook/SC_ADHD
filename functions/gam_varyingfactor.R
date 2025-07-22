# Factor effects will be computed at different smooth_var values.
# int_var should be a factor variable.
gam.varyingcoefficients <- function(region, zscale = F, dataname, smooth_var, int_var, covariates, knots, set_fx = FALSE, increments, draws, return_posterior_coefficients = FALSE){
  
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  #Fit the gam
  gam.data <- get(dataname)
  if (zscale==T){
    gam.data[[region]] <- scale(gam.data[[region]])
  }
  
  smoothmin <- min(gam.data[ , smooth_var])
  smoothmax <- max(gam.data[ , smooth_var])
  
  parcel <- region
  if (is.na(covariates)){
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s)+ %5$s", region, smooth_var, knots, set_fx, int_var))
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s) + %5$s + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  }
  
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame, used to estimate (posterior) model slopes (varying covariate coefficients)
  theseVars <- attr(gam.model$terms,"term.labels") 
  varClasses <- attr(gam.model$terms,"dataClasses") 
  
  #prediction df for int_var level0 at np smooth_var increments
  pred.0 <- data.frame(init = rep(0,np)) 
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == int_var) { 
      #pred.0[,int_var] <- (min(df[,int_var],na.rm = T))
      pred.0[,int_var] <- (levels(df[,int_var])[1])
    } else if (thisVar == smooth_var) {
      pred.0[,smooth_var] = seq(smoothmin,smoothmax, length.out = np)
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {pred.0[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {pred.0[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on first level of factor 
              "ordered" = {pred.0[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on first level of ordinal variable
      )}}
  
  #prediction df for int_var level1 at np smooth_var increments
  pred.1 <- data.frame(init = rep(0,np)) 
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == int_var) { 
      #pred.1[,int_var] <- (max(df[,int_var],na.rm = T))
      pred.1[,int_var] <- (levels(df[,int_var])[2])
    } else if (thisVar == smooth_var) {
      pred.1[,smooth_var] = seq(smoothmin,smoothmax, length.out = np)
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {pred.1[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {pred.1[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on first level of factor 
              "ordered" = {pred.1[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on first level of ordinal variable
      )}}
  
  pred <- rbind(pred.0, pred.1) #complete pred df 
  pred <- pred %>% dplyr::select(-init)
  
  #Get effects (slopes) along the smooth function for the true model
  if(return_posterior_coefficients == FALSE){
    #varying coefficient slopes
    predicted.values <- fitted_values(object = gam.model, data = pred) #predict y at min and mix int_var along the smooth function
    predicted.values <- predicted.values %>% dplyr::select(contains("fitted"), all_of(smooth_var), all_of(int_var))
    colnames(predicted.values) <- c("fitted", "smooth.var", "int.var")
    predicted.values$smooth.var <- round(predicted.values$smooth.var, 3)
    
    varyingcoeff.slopes <-  predicted.values %>% #calculate the effect of int_var on y  (slope; delta y/delta int_var) along the smooth function 
      group_by(smooth.var) %>%
      do(slope = diff(.$fitted)) %>%
      unnest(cols = c(slope))
    colnames(varyingcoeff.slopes) <- c(smooth_var, "slope")  
  }
  
  #Estimate posterior distribution of effects (slopes) from simulated GAM posterior distribution
  if(return_posterior_coefficients == TRUE){
    Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix of fitted gam coefficients
    sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulated model coefficients for npd draws from the posterior
    X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") #get matrix of linear predictors that maps model parameters to the smooth fit (outcome measure scale)
    predicted.values.posterior <- X0 %*% t(sims) #predicted/fitted values along smooth_var, i.e., posterior smooths
    
    predicted.values.posterior <- as.data.frame(predicted.values.posterior)
    colnames(predicted.values.posterior) <- sprintf("draw%s",seq(from = 1, to = npd))
    predicted.values.posterior <- cbind(as.numeric(pred[,smooth_var]), as.factor(pred[,int_var]), predicted.values.posterior)
    colnames(predicted.values.posterior)[1] <- c("smooth.var")
    colnames(predicted.values.posterior)[2] <- c("int.var")
    predicted.values.posterior$smooth.var <- round(predicted.values.posterior$smooth.var, 3)
    
    varyingcoeff.slopes.CI = predicted.values.posterior %>% pivot_longer(cols = contains("draw"), names_to = "draw",values_to = "posterior.fitted") 
    varyingcoeff.slopes.CI$int.var <- factor(varyingcoeff.slopes.CI$int.var, levels=levels(df[,int_var]), labels = c("level0", "level1"))
    varyingcoeff.slopes.CI <- varyingcoeff.slopes.CI %>% pivot_wider(names_from = "int.var", values_from = "posterior.fitted") %>% mutate(slope = (level0 - level1)) 
    #calculate the effect of int_var on y  (slope; delta y/delta int_var) along the smooth function for all draws
    
    varyingcoeff.slopes.CI <- varyingcoeff.slopes.CI %>% dplyr::select(smooth.var, draw, slope, level0, level1)
    varyingcoeff.slopes.CI <- cbind(as.character(parcel), varyingcoeff.slopes.CI)
    colnames(varyingcoeff.slopes.CI) <- c("label", smooth_var, "draw", "slope", levels(df[,int_var])[1], levels(df[,int_var])[2])
  }
  
  if(return_posterior_coefficients == FALSE)
    return(varyingcoeff.slopes)
  if(return_posterior_coefficients == TRUE)
    return(varyingcoeff.slopes.CI)
}

