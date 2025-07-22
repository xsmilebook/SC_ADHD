library(psych)
library(tidyverse)
library(R.matlab)
library(rcompanion)
library(DescTools)
library(polycor)
library(car)
# nominal.vars are nominal variables with more than 2 types.
# dichotomous.vars are nominal variables with 2 types.
corrmat <- function(data, totalvars, dichotomous.vars=NA, ordinal.vars=NA, continuous.vars=NA, nominal.vars=NA){
  corrmat_total <- matrix(NA, nrow=length(totalvars), ncol=length(totalvars))
  colnames(corrmat_total) <- rownames(corrmat_total) <- totalvars
  
  for (x in totalvars){
    for (y in totalvars){
      #print(paste0(x, "-", y))
      var.x <- data[,x]
      var.y <- data[,y]
      test_tab<-table(var.x, var.y)
      rowsum.test <- rowSums(test_tab)
      colsum.test <- colSums(test_tab)
      sum.score <- sum(rowsum.test)
      sum.test <- c(rowsum.test, colsum.test)
      tmp.df <- data.frame(var.x, var.y)
      tmp.df<-tmp.df %>% drop_na()
      
      if (nrow(test_tab)==1 | ncol(test_tab)==1 | sum.score < 1/2*nrow(data) | x==y){
        corr.est <- NA
      }else if (x %in% dichotomous.vars & y %in% dichotomous.vars){
        if (0 %in% sum.test){
          corr.est <- NA
        }else{
          corr.est <- tetrachoric(test_tab)$rho
        }
      }else if ((x %in% dichotomous.vars & y %in% nominal.vars) | (y %in% dichotomous.vars & x %in% nominal.vars) | (y %in% nominal.vars & x %in% nominal.vars)){
        if (0 %in% sum.test){corr.est <- NA}else{
          corr.est <- ContCoef(test_tab, correct = T)
        }
      }else if ((x %in% dichotomous.vars & y %in% ordinal.vars) | (y %in% dichotomous.vars & x %in% ordinal.vars)){
        if (0 %in% sum.test){corr.est <- NA}else{
          if (x %in% dichotomous.vars){
            corr.est <-wilcoxonRG(tmp.df$var.y, tmp.df$var.x)
          }else{
            corr.est <-wilcoxonRG(tmp.df$var.x, tmp.df$var.y)
          }}
      }else if (x %in% ordinal.vars & y %in% ordinal.vars){
        if (0 %in% sum.test){corr.est <- NA}else{
          corr.est <- polycor::polychor(tmp.df$var.x, tmp.df$var.y)}
      }else if (x %in% continuous.vars & y %in% continuous.vars){
        corr.est <- cor(tmp.df, method="pearson")[2,1]
      }else if ((x %in% dichotomous.vars & y %in% continuous.vars) | (y %in% dichotomous.vars & x %in% continuous.vars)){
        if (x %in% dichotomous.vars){
          if (0 %in% rowsum.test){corr.est <- NA}else{
            corr.est <-polycor::polyserial(tmp.df$var.y, tmp.df$var.x)}
        }else{
          if (0 %in% colsum.test){corr.est <- NA}else{
            corr.est <-polycor::polyserial(tmp.df$var.x, tmp.df$var.y)}
        }
      }else if ((x %in% ordinal.vars & y %in% continuous.vars) | (y %in% ordinal.vars & x %in% continuous.vars)){
        if (x %in% ordinal.vars){
          if (0 %in% rowsum.test){corr.est <- NA}else{
            corr.est <-polycor::polyserial(tmp.df$var.y, tmp.df$var.x)}
        }else{
          if (0 %in% colsum.test){corr.est <- NA}else{
            corr.est <-polycor::polyserial(tmp.df$var.x, tmp.df$var.y)}
        }
      }
      corrmat_total[x,y] <- corr.est
    }
  }
  diag(corrmat_total)<-1
  
  return(corrmat_total)
}
