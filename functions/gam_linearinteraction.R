#### Fit GAMM containing interaction term of continuous by categorical variables, which does not contain the smooth var ####
## discrete interaction covariate
gam.linearvar.predict.interaction <- function(region, dataname, smooth_var, VOI, int_var, covariates=NA, knots, set_fx = FALSE){
  gam.data <- get(dataname)
  # tmp<-gam.data[,region]
  # outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  # if (length(outlierindx)>0){
  #   gam.data<-gam.data[-outlierindx, ]
  # }
  parcel <- region
  
  #Fit the gam
  if (is.na(covariates)){
    modelformula <- as.formula(sprintf("%1$s ~ %2$s * %6$s + s(%3$s, k=%4$s, fx=%5$s, by=%6$s)", region, VOI, smooth_var, knots, set_fx, int_var))
    modelformula.null <- as.formula(sprintf("%1$s ~ %2$s + %6$s + s(%3$s, k=%4$s, fx=%5$s, by=%6$s)", region, smooth_var, knots, set_fx, int_var))
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ %2$s * %6$s + s(%3$s, k=%4$s, fx=%5$s, by=%6$s) + %7$s", region, VOI, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null <- as.formula(sprintf("%1$s ~ %2$s + %6$s + s(%3$s, k=%4$s, fx=%5$s, by=%6$s) + %7$s", region, VOI, smooth_var, knots, set_fx, int_var, covariates))
  }
  
  gam.model <- gam(modelformula, method="REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data)
  gam.null.results <- summary(gam.model.null)
  
  # stats
  bootstrap.result.int <- anovaPB(gam.model.null,gam.model, n.sim = 1000,test='Chisq', ncpus=1)
  bootstrap_pvalue <- bootstrap.result.int[1]
  bootstrap_chisq <- bootstrap.result.int[2]
  
  gam.int.pvalue <- gam.results$p.table[grep(x=rownames(gam.results$p.table),pattern = ":"),4]
  # interaction effect size
  sse.model <- sum((gam.model$gam$y - gam.model$gam$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$gam$y - gam.model.null$gam$fitted.values)^2)
  IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  stats.reults<-cbind(parcel, int_var, bootstrap_pvalue, bootstrap_chisq, gam.int.pvalue, IntpartialRsq)
  
  return(stats.reults)
  
}