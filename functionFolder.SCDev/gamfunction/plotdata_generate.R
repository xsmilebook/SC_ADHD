# Modified based on V. J. Sydnor et al. (2023)
# generate 1000 fitted values at 1000 equispaced values of smooth_var for GAM or GAMM models.
# non-interested variables were set as median or mode.
plotdata_generate <- function(modobj, dataname=NA, smooth_var, smoothvector=NA){
  
  if (any(class(modobj)=="gam")) {
    model <- modobj
  } else if (class(modobj$gam)=="gam") {
    model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  
  np <- 1000 #number of predicted values
  if (is.na(dataname)){
    df <- model$model
  }else{
    df <- get(dataname)
  }
  
  theseVars <- attr(model$terms,"term.labels")
  varClasses <- attr(model$terms,"dataClasses")
  thisResp <- as.character(model$terms[[2]])
  # line plot with no interaction
  thisPred <- data.frame(init = rep(0,np))
  
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisClass=="character"){
      df[,thisVar] <- as.factor(df[,thisVar])
      thisClass <- "factor"
    }
    
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
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]}
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init)
  
  #pred$Sex<-levels(Behavior$Sex)[2]
  #pred$handnessfactor<-levels(Behavior$handnessfactor)[3]
  p<-data.frame(predict(model,pred,interval="prediction",level = 0.95, se.fit = T))
  pred <- cbind(pred,p)
  pred$selo <- pred$fit - 1.96*pred$se.fit
  pred$sehi <- pred$fit + 1.96*pred$se.fit
  pred$fit.C <- scale(pred$fit, center = T, scale = F) # zero-centerd fitted values
  pred$fit.Z <- scale(pred$fit, center = T, scale = T) # z-scored fitted values
  pred$fit.floor <- pred$fit-pred$fit[1] # fitted values minus the initial fit
  pred$fit.ratio <- pred$fit / pred$fit[1] # fitted values divided by the initial fit
  pred$selo.C <- pred$fit - 1.96*pred$se.fit-mean(pred$fit)
  pred$sehi.C <- pred$fit + 1.96*pred$se.fit-mean(pred$fit)
  pred[,thisResp] = 1
  return(pred)
}


