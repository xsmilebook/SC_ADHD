library(tidyr)
library(tidyverse)
library(dplyr)


#### FIT GLM logistic ####
##Function to fit a GLM (measure ~ x_var + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
glm.fit.logit <- function(region, dataname, x_var, covariates, mod_only = FALSE){
  
  #Fit the glm
  glm.data <- get(dataname)
  number0 <- sum(glm.data[,region]==0)
  #parcel <- paste0("SC.", as.character(region))
  parcel <- as.character(region)
  modelformula <- as.formula(sprintf("%s ~ %s + %s", region, x_var, covariates))
  glm.model <- glm(modelformula, family=binomial(link="logit"), data = glm.data)
  glm.results <- summary(glm.model)
  
  #glm statistics
  #Z value for the x term and glm-based significance of the x term
  glm.x.Beta <- glm.results$coefficients[2, 1]
  glm.x.se <- glm.results$coefficients[2, 2]
  glm.x.Z <- glm.results$coefficients[2, 3]
  glm.x.pvalue <- glm.results$coefficients[2, 4]
  glm.x.OR <- exp(glm.x.Beta)
  
  stats.results <- cbind(parcel,number0,x_var, glm.x.Beta, glm.x.se, glm.x.Z,glm.x.pvalue, glm.x.OR)
  
  if(mod_only == TRUE){
    return(glm.model)
  }else{
    return(stats.results)
  }
}