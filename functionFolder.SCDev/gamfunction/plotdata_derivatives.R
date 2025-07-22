#### DERIVATIVES ####
##Function to compute smooth derivatives for a main GAM model and for individual draws from the simulated posterior distribution
plot.derivatives <- function(modobj, smooth_var, increments, genderInt=FALSE){
  
  #Set parameters
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  #extract the gam
  if (any(class(modobj)=="gam")) {
    gam.model <- modobj
  } else if (class(modobj$gam)=="gam") {
    gam.model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame, used to estimate (posterior) model coefficients
  thisPred <- data.frame(init = rep(0,np)) 
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to fit for each posterior draw
  for (v in c(1:length(theseVars))) { #fill the prediction df with data 
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init) #prediction df
  
  #Estimate smooth derivatives
  derivs <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = UNCONDITIONAL, data = pred) #derivative at 200 indices of smooth_var with a simultaneous CI
  derivs.fulldf <- derivs %>% dplyr::select(data, derivative, se, lower, upper)
  derivs.fulldf <- derivs.fulldf %>% mutate(significant = !(0 > lower & 0 < upper))
  derivs.fulldf$significant.derivative = derivs.fulldf$derivative*derivs.fulldf$significant
  derivs.fulldf$significant.derivative[derivs.fulldf$significant.derivative==0]<-NA
  colnames(derivs.fulldf) <- c(sprintf("%s", smooth_var), "derivative", "se", "lower", "upper", "significant", "significant.derivative")
  
  return(derivs.fulldf)
}