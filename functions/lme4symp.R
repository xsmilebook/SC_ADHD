library(tidyr)
library(psych)
library(dplyr)
library(lme4)
library(lmerTest)

#### Fit the linear mixed effects of a continuous variable with linear covariates ####
## Function to fit evaluate the linear effects of a continuous variable on dependent variables per edge.
## Interaction effects between continuous variable A and factor or continuous variable B will also be tested.
lmer.fit.symp <- function(region, dataname, symp_var, int_var, covariates=NA, if_residual=F){
  
  #Fit the gam
  lmer.data <- get(dataname)
  symp<-lmer.data[ ,symp_var]
  if (is.numeric(symp)){
    outlierindx1<-which(symp<mean(symp,na.rm = T)-3*sd(symp,na.rm = T) | symp>mean(symp,na.rm = T)+3*sd(symp,na.rm = T))
    lmer.data[outlierindx1 ,symp_var]<-NA
  }
  NonNANIndex <- which(!is.na(lmer.data[ ,symp_var]) & !is.na(lmer.data[ ,region]))
  
  lmer.data <- lmer.data[NonNANIndex,]
  tmp<-lmer.data[,region]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    lmer.data<-lmer.data[-outlierindx, ]
  }
  parcel <- as.character(region)
  
  # Covariates effects
  if (! is.na(covariates)){
    modelformula <- as.formula(sprintf("%s ~ %s + %s + (1|subID)",region, symp_var, covariates))
    modelformula.null<-as.formula(sprintf("%s ~ %s + (1|subID)",region, covariates))
    lmer.model <- lmer(modelformula, data = lmer.data)
    lmer.model.null <- lmer(modelformula.null, data = lmer.data)
    anova.cov.pvalue<-anova(lmer.model, lmer.model.null)$`Pr(>Chisq)`[2]
    Fstat <- anova(lmer.model, lmer.model.null)$Chisq[2]
    #lmer statistics
    lmer.results <- summary(lmer.model)
  }else{
    modelformula <- as.formula(sprintf("%s ~ %s + (1|subID)",region, symp_var))
    lmer.model <- lmer(modelformula, data = lmer.data)
    anova.cov.pvalue<-anova(lmer.model)$`Pr(>F)`[1]
    Fstat <- anova(lmer.model)$`F value`[1]
    #lmer statistics
    lmer.results <- summary(lmer.model)
  }
  
  # Stats in models
  lmer.symp.t <- lmer.results$coefficients[2,4]
  lmer.symp.pvalue <- lmer.results$coefficients[2,5]
  
  # Interaction effects
  if (! is.na(covariates)){
    modelformula.int <- as.formula(sprintf("%1$s ~ %2$s*%3$s + %4$s + (1|subID)",region, symp_var, int_var, covariates))
    modelformula.int.null<-as.formula(sprintf("%1$s ~ %2$s + %3$s + %4$s + (1|subID)",region, symp_var, int_var, covariates))
    lmer.model.int <- lmer(modelformula.int, data = lmer.data)
    lmer.model.int.null <- lmer(modelformula.int.null, data = lmer.data)
    anova.int.pvalue<-anova(lmer.model.int, lmer.model.int.null)$`Pr(>Chisq)`[2]
    Fstat.int <- anova(lmer.model.int, lmer.model.int.null)$Chisq[2]
    #lmer statistics
    lmer.int.results <- summary(lmer.model.int)
  }else{
    modelformula.int <- as.formula(sprintf("%1$s ~ %2$s*%3$s + (1|subID)",region, symp_var, int_var))
    modelformula.int.null<-as.formula(sprintf("%1$s ~ %2$s + %3$s + (1|subID)",region, symp_var, int_var))
    lmer.model.int <- lmer(modelformula.int, data = lmer.data)
    lmer.model.int.null <- lmer(modelformula.int.null, data = lmer.data)
    anova.int.pvalue<-anova(lmer.model.int, lmer.model.int.null)$`Pr(>Chisq)`[2]
    Fstat.int <- anova(lmer.model.int, lmer.model.int.null)$Chisq[2]
    #lmer statistics
    lmer.int.results <- summary(lmer.model.int)
  }
  
  
  #residual correlation
  if (is.na(covariates)){
    res1.normtest <- ks.test(lmer.data[[region]], "pnorm", mean = mean(lmer.data[[region]]), sd = sd(lmer.data[[region]]))
    res2.normtest <- ks.test(lmer.data[[symp_var]], "pnorm", mean = mean(lmer.data[[symp_var]]), sd = sd(lmer.data[[symp_var]]))
    if (res1.normtest$p.value > 0.01 & res2.normtest$p.value > 0.01){
      corrmethod <- "pearson"
    }else{
      corrmethod <- "spearman"
    }
    PCorr_Test <- corr.test(lmer.data[[region]], lmer.data[[symp_var]], method=corrmethod)
    correstimate<-as.numeric(PCorr_Test$r)
    corrp <- as.numeric(PCorr_Test$p)
  }else{
    res1 <- residuals(lmer.model.null)
    modelformula.symp<-as.formula(sprintf("%s ~ %s + (1|subID)",symp_var, covariates))
    lmer.model.symp <- lmer(modelformula.symp, data = lmer.data)
    res2 <- residuals(lmer.model.symp)
    
    res1.normtest <- ks.test(res1, "pnorm", mean = mean(res1), sd = sd(res1))
    res2.normtest <- ks.test(res2, "pnorm", mean = mean(res2), sd = sd(res2))
    if (res1.normtest$p.value > 0.01 & res2.normtest$p.value > 0.01){
      corrmethod <- "pearson"
    }else{
      corrmethod <- "spearman"
    }
    PCorr_Test <- corr.test(res1, res2, method=corrmethod)
    correstimate<-as.numeric(PCorr_Test$r)
    corrp <- as.numeric(PCorr_Test$p)
  }
  
  
  stats.results <- cbind(parcel, symp_var, int_var, lmer.symp.t, lmer.symp.pvalue,anova.cov.pvalue, Fstat, 
                         anova.int.pvalue, Fstat.int, correstimate, corrp)
  residual.sum <- data.frame(SCDeviation.res=res1, Symptom.res=res2, subID=lmer.data$subID)
  
  if (if_residual == F){
    return(stats.results)
  }else{
    return(residual.sum)
  }
  
}

