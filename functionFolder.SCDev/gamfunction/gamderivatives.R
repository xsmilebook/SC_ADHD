#### DERIVATIVES ####
library(gratia)

##Function to compute smooth derivatives for a main GAM model and for individual draws from the simulated posterior distribution
gam.derivatives <- function(modobj,smooth_var, smoothvector=NA, draws, increments, return_posterior_derivatives = TRUE){
  
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior derivative sets estimated
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  EPS <- 1e-07 #finite differences
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  #extract the gam
  if (any(class(modobj)=="gam")) {
    gam.model <- modobj
  } else if (class(modobj$gam)=="gam") {
    gam.model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  region<-gam.model$terms[[2]]
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  theseVars <- attr(gam.model$terms,"term.labels")
  varClasses <- attr(gam.model$terms,"dataClasses")
  thisResp <- as.character(gam.model$terms[[2]])
  #Create a prediction data frame, used to estimate (posterior) model coefficients
  thisPred <- data.frame(init = rep(0,np)) 
  
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) {
      if (is.na(smoothvector[1])){
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      }else{
        thisPred[,smooth_var] = seq(min(smoothvector,na.rm = T),max(smoothvector,na.rm = T), length.out = np)
      }
      
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]},
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]},
              "character" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]}
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init)
  pred2 <- pred #second prediction df
  pred2[,smooth_var] <- pred[,smooth_var] + EPS #finite differences
  
  #Estimate smooth derivatives
  derivs <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = UNCONDITIONAL, data = pred) #derivative at 200 indices of smooth_var with a simultaneous CI
  derivs.fulldf <- derivs %>% dplyr::select(c(smooth_var, ".derivative", ".se", ".lower_ci", ".upper_ci"))
  derivs.fulldf <- derivs.fulldf %>% mutate(significant = !(0 > .lower_ci & 0 < .upper_ci))
  derivs.fulldf$significant.derivative = derivs.fulldf$.derivative*derivs.fulldf$significant
  colnames(derivs.fulldf) <- c(sprintf("%s", smooth_var), "derivative", "se", "lower_ci", "upper", "significant", "significant.derivative")
  
  #Estimate posterior smooth derivatives from simulated GAM posterior distribution
  if(return_posterior_derivatives == TRUE){
    Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
    sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance
    X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") #get matrix of linear predictors for pred
    X1 <- predict(gam.model, newdata = pred2, type = "lpmatrix") #get matrix of linear predictors for pred2
    Xp <- (X1 - X0) / EPS 
    posterior.derivs <- Xp %*% t(sims) #Xp * simulated model coefficients = simulated derivatives. Each column of posterior.derivs contains derivatives for a different draw from the simulated posterior distribution
    posterior.derivs <- as.data.frame(posterior.derivs)
    # posterior.derivs[1,1]<-sum(Xp[1,] * sims[1,]), posterior.derivs[1,2]<-sum(Xp[1,] * sims[2,])...
    # posterior.derivs is the change of posterior smooth term when there is an infinitesimal change of age
    # namely posterior derivatives
    colnames(posterior.derivs) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
    posterior.derivs <- cbind(as.numeric(pred[,smooth_var]), posterior.derivs) #add smooth_var increments from pred df to first column
    colnames(posterior.derivs)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
    posterior.derivs <- cbind(as.character(region), posterior.derivs) #add parcel label to first column
    colnames(posterior.derivs)[1] <- "label_ID" #label the column
    posterior.derivs.long <- posterior.derivs %>% pivot_longer(contains("draw"), names_to = "draw",values_to = "posterior.derivative")
  } #np*npd rows, 3 columns (smooth_var, draw, posterior.derivative)
  
  if(return_posterior_derivatives == FALSE)
    return(derivs.fulldf)
  if(return_posterior_derivatives == TRUE)
    return(posterior.derivs.long)
}