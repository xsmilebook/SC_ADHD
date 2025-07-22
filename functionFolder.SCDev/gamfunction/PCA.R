library(psych)
library(tidyverse)
library(R.matlab)
library(rcompanion)
library(polycor)
library(car)
PCA <- function(data, corrmat_input, scalename, resultFolder, factnum.force=NA){
  ## Step 3 assessing the factorability of the data
  numsam <- sum(complete.cases(data[ ,rownames(corrmat_input)]))
  chisqres<-cortest.bartlett(corrmat_input, n=numsam) # p<0.05
  cortest.bartlett.p<-chisqres$p.value
  KMOres<-KMO(corrmat_input) # Overall MSA >0.49
  overallmsa <- KMOres$MSA
  ## Step 4 screen plot, to determine the number of factors
  #scree(corrmat_input, pc=FALSE)
  faplot<-fa.parallel(corrmat_input, fa="pc", n.obs=numsam)
  if (is.na(factnum.force)){
    factnum<-faplot$ncomp # factor number = 2
  }else{
    factnum<-factnum.force
  }
  pc.out <- principal(corrmat_input, n.obs=numsam, nfactors = 1, max.iter = 100, rotate = "oblimin",
                      covar=T)
  
  ## Step 5 Examine the solution for interpretability
  values.rotation<-pc.out$values[1:factnum] # 10.967587  3.849577
  #100*pc.out$values[1:2]/length(pc.out$values) # 40.62069 14.25769, factor1 explain 40.6% variance, factor2 explain 14.3%
  p <- nrow(pc.out$loadings)
  ssloadings <- colSums(pc.out$loadings ^ 2)
  proportionvar <- colSums(pc.out$loadings ^ 2) / p
  cumulativevar <- cumsum(proportionvar)
  totalexpvar <- cumulativevar[length(cumulativevar)] # the total variationo factors explained
  # ss loadings is the sum of squared loadings. 
  # a factor is worth keeping if the SS loading is greater than 1.
  load <- pc.out$loadings
  write.csv(load, paste0(resultFolder, '/PCA/', scalename, '_loadings', as.character(factnum), 'factor.csv'))
  
  ## calculate cronbach's alpha
  factor.vas <- list()
  alphanum <- rep(NA, times=factnum)
  for (i in 1:factnum){
    factor.vas[[i]]<-rownames(load)[abs(load[,i])>0.3]
    alpha.tmp<-psych::alpha(corrmat_input[factor.vas[[i]],factor.vas[[i]]], check.keys=TRUE)
    alphanum[i]<-alpha.tmp$total$raw_alpha
  }
  
  ## Step 6 calculate factor scores
  data <- data %>% drop_na(colnames(corrmat_input))
  data[,colnames(corrmat_input)] <- lapply(data[,colnames(corrmat_input)], as.numeric)
  factorscores<-factor.scores(data[,colnames(corrmat_input)], pc.out, method="tenBerge")
  factorscores.df <- as.data.frame(factorscores$scores)
  factorscores.df$subID <- data$src_subject_id
  names(factorscores.df)[1:factnum]<-paste0(scalename, c(1:factnum))
  
  
  FAinfo <- data.frame(scalename, numsam, cortest.bartlett.p,
                       overallmsa, factnum, totalexpvar)
  FAperformance <- data.frame(ssloadings, proportionvar, cumulativevar)
  FAperformance$factor <- c(1:factnum)
  FAperformance$alpha <- alphanum
  FAtotal.result <- list(pc.out=pc.out, values.rotation=values.rotation, 
                         FAperformance=FAperformance,FAinfo=FAinfo, load=load, factor.vas=factor.vas, factorscores.df=factorscores.df)
  
  return(FAtotal.result)
}