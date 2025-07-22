library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(lme4)
library(gamm4)
library(pbkrtest)
#library(car)

pbootint <- function(modelobj, int_var=NA){
  numsims <- 1000
  set.seed(925)
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- attr(terms(f1),"term.labels")
  
  if (sum(str_detect(theseVars, ":"))==1){
    intindx <- which(str_detect(theseVars, ":"))
    f2 <- reformulate(theseVars[-intindx],response = thisResp)
  }else{
    f2 <- reformulate(theseVars[2:(length(theseVars))],response = thisResp)
  }
  
  g1 <- gam(f1,data = df)
  g2 <- gam(f2,data = df)
  
  mat1 <- model.matrix(g1)
  mat2 <- model.matrix(g2)
  
  subID<-df$subID
  y <- df[,thisResp]
  
  m1 <- lmer(y ~ -1 + mat1 + (1|subID))
  m2 <- lmer(y ~ -1 + mat2 + (1|subID))
  refdist <- PBrefdist(m1, m2, nsim=numsims)
  pb <- PBmodcomp(m1, m2, ref = refdist, details = 1)
  int_pval <- pb$test["PBtest","p.value"]
  int_chisq <- pb$test["PBtest","stat"]
  
  return(c(int_pval, int_chisq))
}

#### Fit GAMM containing interaction term of continuous by categorical variables, which does not contain the smooth var ####
## discrete interaction covariate
gamm.linearvar.predict.interaction <- function(region, dataname, smooth_var, VOI, int_var, covariates=NA, knots, set_fx = FALSE){
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
  
  gamm.model <- gamm4(modelformula, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.results <- summary(gamm.model$gam)
  
  gamm.model.null <- gamm4(modelformula.null, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.null.results <- summary(gamm.model.null$gam)
  
  # stats
  bootstrap.result.int <- pbootint(gamm.model, int_var)
  bootstrap_pvalue <- bootstrap.result.int[1]
  bootstrap_chisq <- bootstrap.result.int[2]
  
  gam.int.pvalue <- gamm.results$p.table[grep(x=rownames(gamm.results$p.table),pattern = ":"),4]
  # interaction effect size
  sse.model <- sum((gamm.model$gam$y - gamm.model$gam$fitted.values)^2)
  sse.nullmodel <- sum((gamm.model.null$gam$y - gamm.model.null$gam$fitted.values)^2)
  IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  stats.reults<-cbind(parcel, int_var, bootstrap_pvalue, bootstrap_chisq, gam.int.pvalue, IntpartialRsq)
  
  return(stats.reults)
  
}
