library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(lme4)
library(gamm4)
library(pbkrtest)

pbootint <- function(modelobj, int_var=NA){
  numsims <- 1000
  set.seed(925)
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- attr(terms(f1),"term.labels")
  if (!is.na(int_var)){
    indexcomma <- gregexpr(",", theseVars[1])[[1]]
    addsmooth <- theseVars[1]
    addsmooth <- paste0(substr(addsmooth, 1, indexcomma[1]), substr(addsmooth, indexcomma[2]+1, nchar(addsmooth)))
    theseVars <-  c(theseVars, addsmooth)
    }
  if (sum(str_detect(theseVars, "by ="))==1){
    int_var <- str_split_i(theseVars[1], "by = ", 2)
    int_var <- str_split_i(int_var, ", ", 1)
    f2 <- reformulate(c(theseVars[2:(length(theseVars))], int_var),response = thisResp)
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

#### PREDICT GAM SMOOTH FITTED VALUES FOR A SPECIFIED VALUE OF AN INTERACTING COVARIATE ####
## discrete interaction covariate
##Function to predict fitted values of a region for a given value of a covariate, using a varying coefficients smooth-by-linear covariate interaction
gamm.smooth.predict.interaction <- function(region, dataname, smooth_var, int_var, covariates=NA, knots, set_fx = FALSE){
  gam.data <- get(dataname)
  # tmp<-gam.data[,region]
  # outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  # if (length(outlierindx)>0){
  #   gam.data<-gam.data[-outlierindx, ]
  # }
  parcel <- region
 
  #Fit the gam
  if (is.na(covariates)){
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s)+ %5$s", region, smooth_var, knots, set_fx, int_var))
    modelformula.null <- as.formula(sprintf("%1$s ~ %5$s+s(%2$s, k=%3$s, fx=%4$s)", region, smooth_var, knots, set_fx, int_var))
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s) + %5$s + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
    modelformula.null <- as.formula(sprintf("%1$s ~ %5$s+s(%2$s, k=%3$s, fx=%4$s)+ %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  }
  gamm.model <- gamm4(modelformula, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.results <- summary(gamm.model$gam)
  
  gamm.model.null <- gamm4(modelformula.null, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.null.results <- summary(gamm.model.null$gam)
  
  # stats
  T.disease <- gamm.null.results$p.table[2:nlevels(gam.data[[int_var]]),3]
  P.disease <- gamm.null.results$p.table[2:nlevels(gam.data[[int_var]]),4]
  bootstrap.result.disease <- pbootint(gamm.model.null)
  bootstrap.P.disease <- bootstrap.result.disease[1]
  bootstrap.chisq.disease <- bootstrap.result.disease[2]
  #bootstrap.P.disease<-NA
  
  bootstrap.result.int <- pbootint(gamm.model, int_var)
  bootstrap_pvalue <- bootstrap.result.int[1]
  bootstrap_chisq <- bootstrap.result.int[2]
  
  gam.int.pvalue <- gamm.results$s.table[grep(x=rownames(gamm.results$s.table),pattern = int_var),"p-value"]
  # interaction effect size
  sse.model <- sum((gamm.model$gam$y - gamm.model$gam$fitted.values)^2)
  sse.nullmodel <- sum((gamm.model.null$gam$y - gamm.model.null$gam$fitted.values)^2)
  IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  T.disease <- c(1, T.disease)
  P.disease <- c(NA, P.disease)
  stats.reults<-cbind(parcel, int_var, bootstrap_pvalue, bootstrap_chisq, gam.int.pvalue, IntpartialRsq, T.disease, 
                      P.disease, bootstrap.P.disease, bootstrap.chisq.disease)
  
  return(stats.reults)

}






