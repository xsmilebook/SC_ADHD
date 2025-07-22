library(tidyverse)
library(gamlss)
library(scales)
library(parallel)

### Function 1
# This function is used to compare GAMLSS with different distribution.
# Referring to Bethlehem et. al., Nature, 2022 & Sun et. al., bioArxiv, 2023,
# we firstly compared the fitting performance of different distributions to select the best one.
# Only continuous distributions with >= 3 parameters were compared.

gamlss_comparedistribution <- function(dataname, dependentvar, smoothvar,IDvar, bs.df, covariates, randomvar=NA, cl){
  
  # get data
  gam.data <- get(dataname)
  gam.data$age <- gam.data[[smoothvar]]
  
  covariates <- gsub(" ", "", covariates)
  if (! is.na(randomvar)){
    gam.data1 <- gam.data %>% select(c(dependentvar, "age", unlist(strsplit(covariates, "+", fixed=T)), randomvar, IDvar)) %>% drop_na()
    gam.data1$dependent <- gam.data1[[dependentvar]]
    gam.data1$randomvar <- as.factor(gam.data1[[randomvar]])
    #gam.data1 <- gam.data1 %>% filter((dependent>mean(dependent)-3*sd(dependent)) & (dependent<mean(dependent)+3*sd(dependent)))
    gam.data1$dependent <- NULL
    gam.data2 <- gam.data1 %>% group_by(randomvar) %>%
      filter(n() > 30) %>%
      ungroup()
    
    mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s) + %s + random(as.factor(%s))", dependentvar, "age", bs.df, covariates, "randomvar"))
  }else{
    gam.data2 <- gam.data %>% select(c(dependentvar, "age", unlist(strsplit(covariates, "+", fixed=T)), IDvar)) %>% drop_na()
    mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s) + %s", dependentvar, "age", bs.df, covariates))
  }
  
  # fit models. 26 distributions were adopted.
  #	names	No_parameters	continuous
  #	BEINF()	4	1
  #	GG()	3	1
  #	GIG()	3	1
  #	GT()	4	1
  #	PE()	3	1
  #	PE2()	3	1
  #	SEP1()	4	1
  #	SEP2()	4	1
  #	SEP3()	4	1
  #	SEP4()	4	1
  #	SHASH()	4	1
  #	JSU()	4	1
  #	SHASHo()	4	1
  #	ST1()	3	1
  #	ST2()	3	1
  #	ST3()	3	1
  #	ST4()	3	1
  #	ST5()	3	1
  #	TF()	3	1
  #	BCCG()	3	1
  #	BCPE()	4	1
  #	BCT()	4	1
  #	exGAUS()	3	1
  #	EGB2()	4	1
  #	GB1()	4	1
  #	GB2()	4	1
  
  
  familylist <- c("JSU", "BCCG", "BCPE", "BCT", "exGAUS", "EGB2", 
                  "GB2", "GG", "GIG", "GT", "PE", "PE2", "SEP1", "SEP2", 
                  "SEP3", "SEP4", "SHASH", "SHASHo", "ST1", "ST2", "ST3", 
                  "ST4", "ST5", "TF")
  con<-gamlss.control(n.cyc=200)
  
  clusterExport(cl, varlist = ls(), envir = environment())
  mod.sum <- parLapply(cl, 1:length(familylist), function(i) {
  
    result <- try({
      command <- paste0(
        "mod.tmp <- gamlss(mod.mu.formula, sigma.formula =~ bs(age) + ",
        covariates, 
        ", nu.formula = ~1, family=", familylist[i], ", data=gam.data2, control=con)"
       )
    
      eval(parse(text = command))
    
      mod.tmp$ID <- gam.data2[[IDvar]]
      return(mod.tmp)
   }, silent = TRUE)
  
     if (inherits(result, "try-error")) {
       return(list(error = TRUE, message = attr(result, "condition")$message))
     }
  
     return(result)
  })


  # mod.sum <- list()
  # for (i in 1:length(familylist)){
  #   command <- paste0("mod.sum[[i]] <- gamlss(mod.mu.formula, sigma.formula =~ bs(age) + ", 
  #                     covariates, 
  #                ", nu.formula = ~1,family=", familylist[i],", data=gam.data2, control=con)")
  #   
  #   eval(parse(text = command))
  # }
  
  performance <- data.frame(matrix(NA, length(familylist), 4))
  names(performance) <- c("model_ID", "distribution", "converged", "BIC")
    
  for (i in 1:length(mod.sum)){
    if (length(mod.sum[[i]])>30){
      performance$model_ID[i] <- i
      performance$distribution[i] <- mod.sum[[i]]$family[[1]]
      performance$converged[i] <- mod.sum[[i]]$converged
      performance$BIC[i] <- mod.sum[[i]]$sbc
    }
  }
  familylist <- familylist[performance$converged==T]
  print(paste0("The best distibution for ", dependentvar, " is ", familylist[which.min(performance$BIC[performance$converged==T])], "."))
  
  performance$dependentvar <- dependentvar
  
  outlist <- list(mod.sum, performance=performance)
  
  return(outlist)
  
}


### Function 2
# This function is used to compare GAMLSS with different bs degree of freedom.
# Referring to Sun et. al., bioArxiv, 2023, df = 3~6 for mu & sigma were tested.
# bs.df.mat should be a n*2 matrix with the first column is the df of mu and the second column 
# is the df of sigma.

gamlss_compare.bs.df <- function(dataname, dependentvar, smoothvar, IDvar, bs.df.set, covariates, distribution.fam, randomvar=NA, cl){
  
  # get data
  gam.data <- get(dataname)
  gam.data$age <- gam.data[[smoothvar]]
  covariates <- gsub(" ", "", covariates)
  if (! is.na(randomvar)){
    gam.data1 <- gam.data %>% select(c(dependentvar, "age", unlist(strsplit(covariates, "+", fixed=T)), randomvar, IDvar)) %>% drop_na()
    gam.data1$dependent <- gam.data1[[dependentvar]]
    gam.data1$randomvar <- as.factor(gam.data1[[randomvar]])
    #gam.data1 <- gam.data1 %>% filter((dependent>mean(dependent)-3*sd(dependent)) & (dependent<mean(dependent)+3*sd(dependent)))
    gam.data1$dependent <- NULL
    gam.data2 <- gam.data1 %>% group_by(randomvar) %>%
      filter(n() > 30) %>%
      ungroup()
  }else{
    gam.data2 <- gam.data %>% select(c(dependentvar, "age", unlist(strsplit(covariates, "+", fixed=T)))) %>% drop_na()
  }
  
  con<-gamlss.control(n.cyc=200)
  clusterExport(cl, varlist = ls(), envir = environment())
  mod.sum <- parLapply(cl, 1:nrow(bs.df.set), function(i){
    library(gamlss)
    # df for mu // df for sigma
    mu.df <- bs.df.set[i, 1]
    sigma.df <- bs.df.set[i, 2]
    degree <- bs.df.set[i, 3]
    
    if (! is.na(randomvar)){
      mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree = %s) + %s + random(as.factor(%s))", dependentvar, "age", mu.df, degree, covariates, "randomvar"))
    }else{
      mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree = %s) + %s", dependentvar, "age", mu.df, degree, covariates))
    }
    
    result <- try({
      command <- paste0("mod.tmp <- gamlss(mod.mu.formula, sigma.formula =~ bs(age, df = ", sigma.df, ", degree = ", degree,") + ", 
                        covariates, 
                        ", nu.formula = ~1, family=", distribution.fam,", data=gam.data2, control=con)")
      
      eval(parse(text = command))
      mod.tmp$ID <- gam.data2[[IDvar]]
      return(mod.tmp)
    }, silent = TRUE)
    
    if(inherits(result, "try-error")){
      return(list(error = TRUE, message = attr(result, "condition")$message))
    }
    
    return(mod.tmp)
  })
  
  performance <- data.frame(matrix(NA, nrow(bs.df.set), 7))
  names(performance) <- c("model_ID", "distribution", "converged", "BIC", "mu.df", "sigma.df", "degree")
  
  for (i in 1:nrow(bs.df.set)){
    model_result <- mod.sum[[i]]
    
    # judge whether the model is converged
    if (!is.null(model_result$family)) { # 一个更稳健的检查是看 'family' 元素是否存在
      performance$model_ID[i] <- i
      performance$distribution[i] <- model_result$family[[1]]
      performance$converged[i] <- model_result$converged
      performance$BIC[i] <- model_result$sbc
      performance$mu.df[i] <- bs.df.set[i,1]
      performance$sigma.df[i] <- bs.df.set[i,2]
      performance$degree[i] <- bs.df.set[i,3]
    } else {
      # failed model
      performance$model_ID[i] <- i
      performance$distribution[i] <- distribution.fam # 仍然记录是为哪个分布尝试的
      performance$converged[i] <- FALSE # 明确标记为未收敛
      performance$BIC[i] <- NA # BIC不可用
      performance$mu.df[i] <- bs.df.set[i,1]
      performance$sigma.df[i] <- bs.df.set[i,2]
      performance$degree[i] <- bs.df.set[i,3]
    }
  }
  # for (i in 1:nrow(bs.df.set)){
  #   performance$model_ID[i] <- i
  #   performance$distribution[i] <- mod.sum[[i]]$family[[1]]
  #   performance$converged[i] <- mod.sum[[i]]$converged
  #   performance$BIC[i] <- mod.sum[[i]]$sbc
  #   performance$mu.df[i] <- bs.df.set[i,1]
  #   performance$sigma.df[i] <- bs.df.set[i,2]
  #   performance$degree[i] <- bs.df.set[i,3]
  # }
  
  # bs.df.set <- bs.df.set[performance$converged==T, ]
  print(paste0("The best distibution for ", dependentvar, " is mu.df = ", bs.df.set[which.min(performance$BIC[performance$converged==T]), 1], 
               "sigma.df = ", bs.df.set[which.min(performance$BIC[performance$converged==T]), 2], ", degree = ", bs.df.set[which.min(performance$BIC[performance$converged==T]), 3], "."))
  #saveRDS(mod.sum, paste0(saveout_dir, "/", dependentvar, "_gamlss_26distributions_bsds_", bs.df, ".rds"))
  
  performance$dependentvar <- dependentvar
  outlist <- list(mod.sum, performance=performance)
  
  return(outlist)

}


Boot.Function <- function( n, Base.Seed, dataname,smoothvar,randomvar, model_obj, stratify=NULL ) {
  # extract variables
  mod.mu.formula <- model_obj$mu.formula
  formula.vars <- as.character(mod.mu.formula)
  dependentvar <- formula.vars[2]
  
  covariates <- str_split(formula.vars[3], "\\+")
  if (! is.na(randomvar)){
    covariates <- covariates[[1]][2:(length(covariates[[1]])-1)]
  }else{
    covariates <- covariates[[1]][2:length(covariates[[1]])]
  }
  covariates <- gsub(" ", "", covariates)
  # get data
  gam.data <- get(dataname)
  
  if (! is.na(randomvar)){
    gam.data1 <- gam.data %>% select(unique(c(dependentvar, smoothvar, unlist(strsplit(covariates, "+", fixed=T)), randomvar, stratify))) %>% drop_na()
    gam.data1$dependent <- gam.data1[[dependentvar]]
    gam.data1$randomvar <- as.factor(gam.data1[[randomvar]])
    #gam.data1 <- gam.data1 %>% filter((dependent>mean(dependent)-3*sd(dependent)) & (dependent<mean(dependent)+3*sd(dependent)))
    gam.data1$dependent <- NULL
    gam.data2 <- gam.data1 %>% group_by(randomvar) %>%
      filter(n() > 30) %>%
      ungroup()
    
    gam.data2[[randomvar]] <- droplevels(gam.data2[[randomvar]])
  }else{
    gam.data2 <- gam.data %>% select(c(dependentvar, "age", covariates)) %>% drop_na()
  }
  
  set.seed( seed=Base.Seed + n )
  if (length(stratify) == 1) {
    stratify.var <- gam.data2[[stratify]]
  } else{
    stratify.var <- interaction(gam.data2[ ,stratify])
  }
  
  
  if (length(stratify)==0){
    INDEX <- sample(1:NROW(gam.data2),NROW(gam.data2),replace=TRUE) ## unstratified bootstrap
  } else {
    
    
    SPLIT <- split(1:NROW(gam.data2), stratify.var)
    LAPPLY <- lapply(SPLIT,function(X){sample(x=X,size=length(X),replace=TRUE)})
    INDEX <- unsplit(LAPPLY, stratify.var)

  }
  
  gam.data.subset <- gam.data2[INDEX, ] ## generate SUBSET with bootstrapped-SUBSET

  # fit gamlss
  sigma.formula <- as.character(model_obj$sigma.formula)[2]
  mod.mu.formula <- paste(formula.vars[2], "~", formula.vars[3])
  distribution.fam <- model_obj$family[1]
  con<-gamlss.control(n.cyc=200)
  command <- paste0("mod.tmp <- gamlss(", mod.mu.formula, ", sigma.formula =~", sigma.formula,
                    ", nu.formula = ~1,family=", distribution.fam,", data=gam.data.subset, control=con)")
  #print(summary(gam.data.subset))
  eval(parse(text = command))
  bootstrap_list <- list(gam.data.subset=gam.data.subset, mod.tmp=mod.tmp)
  
  return(bootstrap_list)
}







