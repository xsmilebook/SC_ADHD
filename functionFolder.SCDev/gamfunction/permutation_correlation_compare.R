# This function is used to perform permutation test on the 
# group-wise difference of alignment with connectional axis.
library(tidyverse)
library(lme4)
library(gamm4)
library(mgcv)
library(gratia)
coefficientdif_calc <- function(dataset, regionlist, SCrank.df){
  # gam model to extract real correlation
  gamresult <- data.frame(parcel=regionlist, partialRsq.td=rep(NA, length(regionlist)), 
                          partialRsq.disease=rep(NA, length(regionlist)))
  gamresult <- gamresult %>% left_join(SCrank.df, by="parcel")
  print("1")
  for (region in regionlist){
    if (str_detect(disease, "pFactor")==FALSE){
      index.disease <- which(dataset[,disease]==1)
      index.td <- which(dataset[,"td"]==1)
    }else{
      index.disease <- which(dataset[,disease]==1)
      index.td <- which(dataset[,disease]==0)
    }
    
    gam.data.disease <- dataset[index.disease,]
    gam.data.td <- dataset[index.td,]
    if (interest_var==smooth_var){
      modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
      theseVars <- attr(terms(modelformula),"term.labels")
      modelformula.null <- reformulate(theseVars[2:(length(theseVars))],response = region)
      # td
      gamm.model.td <- gamm4(modelformula, random=~(1|subID), REML=TRUE, data = gam.data.td)
      gamm.model.td.null <- gamm4(modelformula.null, random=~(1|subID), REML=TRUE, data = gam.data.td)
      sse.model <- sum((gamm.model.td$gam$y - gamm.model.td$gam$fitted.values)^2)
      sse.nullmodel <- sum((gamm.model.td.null$gam$y - gamm.model.td.null$gam$fitted.values)^2)
      partialRsq.td <- (sse.nullmodel - sse.model)/sse.nullmodel
      derv <- derivatives(gamm.model.td$gam, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = F)
      mean.derivative <- mean(derv$derivative)
      if (mean.derivative < 0){
        gamresult$partialRsq.td[gamresult$parcel==region] = partialRsq.td*-1
      }else{gamresult$partialRsq.td[gamresult$parcel==region] = partialRsq.td}
      
      # disease
      gamm.model.disease <- gamm4(modelformula, random=~(1|subID), REML=TRUE, data = gam.data.disease)
      gamm.model.disease.null <- gamm4(modelformula.null, random=~(1|subID), REML=TRUE, data = gam.data.disease)
      sse.model <- sum((gamm.model.disease$gam$y - gamm.model.disease$gam$fitted.values)^2)
      sse.nullmodel <- sum((gamm.model.disease.null$gam$y - gamm.model.disease.null$gam$fitted.values)^2)
      partialRsq.disease <- (sse.nullmodel - sse.model)/sse.nullmodel
      derv <- derivatives(gamm.model.disease$gam, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = F)
      mean.derivative <- mean(derv$derivative)
      if (mean.derivative < 0){
        gamresult$partialRsq.disease[gamresult$parcel==region] = partialRsq.disease*-1
      }else{gamresult$partialRsq.disease[gamresult$parcel==region] = partialRsq.disease}
      
    }else{
      # clean data
      cognition<-gam.data.td[ ,interest_var]
      outlierindx1<-which(cognition<mean(cognition,na.rm = T)-3*sd(cognition,na.rm = T) | cognition>mean(cognition,na.rm = T)+3*sd(cognition,na.rm = T))
      gam.data.td[outlierindx1 ,interest_var]<-NA
      NonNANIndex <- which(!is.na(gam.data.td[ ,interest_var]) & !is.na(gam.data.td[ ,region]))

      gam.data.td <- gam.data.td[NonNANIndex,]
      tmp<-gam.data.td[,region]
      outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
      if (length(outlierindx)>0){
        gam.data.td<-gam.data.td[-outlierindx, ]
      }

      cognition<-gam.data.disease[ ,interest_var]
      outlierindx1<-which(cognition<mean(cognition,na.rm = T)-3*sd(cognition,na.rm = T) | cognition>mean(cognition,na.rm = T)+3*sd(cognition,na.rm = T))
      gam.data.disease[outlierindx1 ,interest_var]<-NA
      NonNANIndex <- which(!is.na(gam.data.disease[ ,interest_var]) & !is.na(gam.data.disease[ ,region]))

      gam.data.disease <- gam.data.disease[NonNANIndex,]
      tmp<-gam.data.disease[,region]
      outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
      if (length(outlierindx)>0){
        gam.data.disease<-gam.data.disease[-outlierindx, ]
      }
      # fit model
      modelformula <- as.formula(sprintf("%s ~ %s+s(%s, k= %s, fx = %s) + %s",interest_var, region, smooth_var, knots, set_fx,
      covariates)) 
      theseVars <- attr(terms(modelformula),"term.labels") 
      modelformula.null <- reformulate(theseVars[2:(length(theseVars))],response = interest_var)
      # td
      gam.model.td <- gam(modelformula, method = "REML", data = gam.data.td)
      gam.result.td <- summary(gam.model.td)
      gamresult$partialRsq.td[gamresult$parcel==region] = gam.result.td$p.table[2,3]
      # disease
      gam.model.disease <- gam(modelformula, method = "REML", data = gam.data.disease)
      gam.result.disease <- summary(gam.model.disease)
      gamresult$partialRsq.disease[gamresult$parcel==region] = gam.result.disease$p.table[2,3]
    }}
  # real correlation difference
  coefficient_td <- cor(gamresult$partialRsq.td, gamresult$SCrank, method = "spearman")
  coefficient_disease <- cor(gamresult$partialRsq.disease, gamresult$SCrank, method = "spearman")
  coefficient_dif <- coefficient_disease-coefficient_td
  
  return(c(coefficient_td, coefficient_disease, coefficient_dif))}

perm.alignwithFixedrank <- function(dataname,disease, smooth_var, interest_var, covariates, knots,ds.resolution=12,permnum=1000, set_fx = FALSE, stats_only=TRUE){
  
  #connectional rank
  Matrix.ds<-matrix(NA, nrow=ds.resolution, ncol=ds.resolution)
  indexsave.ds <- lower.tri(Matrix.ds, diag=T)
  Matrix.ds.index<-Matrix.ds
  SClength.ds=ds.resolution*(ds.resolution+1)/2
  Matrix.ds.index[indexsave.ds]<-c(1:SClength.ds)
  Matrix.ds.SCrank<-Matrix.ds
  for (x in 1:ds.resolution){
    for (y in 1:ds.resolution){
      Matrix.ds.SCrank[x,y]<-x^2+y^2
    }
  }
  SCrank<-rank(Matrix.ds.SCrank[indexsave.ds], ties.method = "average")
  SCrank.df <- data.frame(SCrank=SCrank, parcel=paste0("SC.", c(1:SClength.ds), "_h"))
  
  # function calculate coefficient difference
  gam.data <- get(dataname)
  regionlist <- names(gam.data)[str_detect(names(gam.data), "SC.")]
  if (str_detect(disease, "pFactor")==FALSE){
    index.disease <- which(gam.data[,disease]==1 & gam.data[,paste0(disease,"_consistent")]==1)
    index.td <- which(gam.data[,"td"]==1 & gam.data[,"td_consistent"]==1)
  }else{
    index.disease <- which(gam.data[,disease]==1)
    index.td <- which(gam.data[,disease]==0)
  }
  gam.data.disease <- gam.data[index.disease,]
  gam.data.td <- gam.data[index.td,]

  gam.data.merge <- rbind(gam.data.td, gam.data.disease)
  
  ## real coefficient difference
  coefficient.real <- coefficientdif_calc(gam.data.merge, regionlist, SCrank.df)
  coefficient.real <- data.frame(coef.td=coefficient.real[1], coef.disease=coefficient.real[2],
                                    coef.diff=coefficient.real[3])
  coefficient_dif.real <- coefficient.real$coef.diff[1]
  print(coefficient_dif.real)
  ## permutate diagnose
  set.seed(NULL)
  coeficient_dif.null <- mclapply(1:permnum, function(x){
    perm_index <- sample(1:nrow(gam.data.merge))
    gam.data.merge[,disease] <- gam.data.merge[perm_index,disease]
    coefficient_dif.null.tmp <- coefficientdif_calc(gam.data.merge, regionlist, SCrank.df)
    return(coefficient_dif.null.tmp)
  }, mc.cores = 50)
  coeficient_dif.null.vector <- lapply(coeficient_dif.null, function(x) unlist(x)[3])
  coeficient_dif.null.vector <- as.numeric(unlist(coeficient_dif.null.vector))
  print(coeficient_dif.null)
  print(coeficient_dif.null.vector)
  if (coefficient_dif.real>0){
    p_value <- sum(coeficient_dif.null.vector > coefficient_dif.real) / permnum
  }else{p_value <- sum(coeficient_dif.null.vector < coefficient_dif.real) / permnum}
  
  coefficient.real$p_value = p_value
  coefficient.real$interest_var = interest_var
  coefficient.real$disease = disease
  if (stats_only==TRUE){
    return(coefficient.real)
  }else{
    return(list(coefficient.real, coeficient_dif.null.vector))
  }
  
}



