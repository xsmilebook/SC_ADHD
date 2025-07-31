rm(list=ls())
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(MatchIt)
# set resolution
Yeoresolution <- 200

element_num <- 252
# input directory
homepath <- "D:/code/SC_ADHD"

interfileFolder <- file.path(homepath, "datasets", "interfileFolder", "TD_ADHD_deviation")
functionFolder <- file.path(homepath, "code", "functions")
resultFolder <- file.path(homepath, "datasets", "results", "TD_ADHD_deviation")
functionFolder.SCDev <- file.path(homepath, "code", "functionFolder.SCDev", "gamfunction")
FigureFolder <- file.path(homepath, "Figures", "TD_ADHD_deviation")


# 1. read data
SC_labels <- readLines(file.path(interfileFolder, "roi.txt"), ok = TRUE, skipNul = FALSE)

SCdata.train <- read.csv(file.path(interfileFolder, "df_tr.csv"))
SCdata.test <- read.csv(file.path(interfileFolder, "df_te.csv"))


# data clean
SCdata.train$sex <- factor(SCdata.train$sex)
SCdata.test$sex <- factor(SCdata.test$sex)

all_site_levels <- unique(c(as.character(SCdata.train$site), as.character(SCdata.test$site)))

SCdata.train$site <- factor(SCdata.train$site, levels = all_site_levels)
SCdata.test$site <- factor(SCdata.test$site, levels = all_site_levels)



# source function
source(paste0(functionFolder, "/Construct_gamlss_set.R"))


## Fit normative models
dataname <- "SCdata.train"
smoothterm <- "age"
covariates <- "sex+mean_fd+TBV+network_strength_schaefer200"
randomvar <- "site"
mu.df <- sigma.df <- degree <- 2
distribution.fam <- "GG"
IDvar <- "sub_id"
quantile.vec <- c(0.025, 0.5, 0.975)
stratify <- c("sex", "site")
# MARK: debug


if (! file.exists(paste0(interfileFolder, "/GAMLSS_Schaefer200.TDtraining.sum.rds"))){
  cl <- makeCluster(24)
  clusterExport(cl, varlist = ls(), envir = environment())
  invisible(clusterEvalQ(cl, {
    library(ggplot2)
    library(tidyverse)
    library(openxlsx)
    library(parallel)
    library(gamlss)
    library(scales)
    source(paste0(functionFolder, "/Construct_gamlss_set.R"))
  }))
  mod78.training.sum <- parLapply(cl, 1:element_num, function(i){
    tryCatch({
      dependentvar <- SC_labels[i]
      sumlist <- construct_gamlss(dataname, dependentvar, smoothterm, covariates, randomvar, 
                                  mu.df, sigma.df, degree, distribution.fam, IDvar, quantile.vec, stratify)
      
      return(sumlist)
    }, error = function(e){

      error_info <- list(
        status = "error",
        index = i,
        dependent_variable = SC_labels[i],
        error_message = e$message
      )
      
      
      return(error_info)
    })
  })
  saveRDS(mod78.training.sum, paste0(interfileFolder, "/GAMLSS_Schaefer200.TDtraining.sum.rds"))
}else{
  mod78.training.sum <- readRDS(paste0(interfileFolder, "/GAMLSS_Schaefer200.TDtraining.sum.rds"))
}


# judge converged
converged_flag <- TRUE

for(i in 1:element_num){
  performance <- mod78.training.sum[[i]]$performance.tb
  if(is.null(performance)){
    cat(paste0(SC_labels[i], "is not converged! The index is: ", i))
    converged_failure_idx <- i
    converged_flag <- FALSE
    
    # rerun gamlss
    dependentvar <- SC_labels[converged_failure_idx]
    # train a simple model
    mod.mu.simple <- as.formula(paste0(dependentvar, " ~ ", covariates))
    sigma.simple <- as.formula("~ 1")
    simple_family <- GG
    mod_simple <- try({
      gamlss(mod.mu.simple,
             sigma.formula = sigma.simple,
             nu.formula = ~1,
             family = simple_family,
             data = SCdata.train,
             control = gamlss.control(n.cyc = 100, t.tol = 0.001))
    })
    
    
    sumlist <- construct_gamlss(dataname, dependentvar, smoothterm, covariates, randomvar, 
                                mu.df, sigma.df, degree, distribution.fam, IDvar, quantile.vec, stratify, mod_simple)
    mod78.training.sum[[i]] <- sumlist
  }
}
if(!converged_flag){
  saveRDS(mod78.training.sum, paste0(interfileFolder, "/GAMLSS_Schaefer200.TDtraining.sum.rds"))
}



# compute deviations
 
SCdata.requireZ <- SCdata.test

# compute deviation & draw plots
SCdata.requireZ_fixed <- SCdata.requireZ
SCdata.requireZ_fixed$mean_fd <- mean(SCdata.requireZ_fixed$mean_fd)
SCdata.requireZ_fixed$TBV <- mean(SCdata.requireZ_fixed$TBV)
SCdata.requireZ_fixed$network_strength_schaefer200 <- mean(SCdata.requireZ_fixed$network_strength_schaefer200)
SCdata.requireZ_fixed.F <- SCdata.requireZ_fixed.M <- SCdata.requireZ_fixed
SCdata.requireZ_fixed.F$sex <- factor(0)
SCdata.requireZ_fixed.M$sex <- factor(1)

sitelist <- unique(SCdata.train$site)



# fix the bug of gamlss.predict: we must define the data used for training(the dataname gam.data2 is used in the construct_gamlss function) 
gam.data2 <- SCdata.train


if (! file.exists(paste0(interfileFolder, "/SCdata.Schaefer200.testset_ADHD_deviation.rds"))){
  cl <- makeCluster(24)
  
  clusterExport(cl, varlist = ls(), envir = environment())
  invisible(clusterEvalQ(cl, {
    library(ggplot2)
    library(tidyverse)
    library(openxlsx)
    library(parallel)
    library(gamlss)
    library(scales)
    source(paste0(functionFolder, "/Construct_gamlss_set.R"))
  }))
  clusterEvalQ(cl, {
    sink("worker_output.log", append = TRUE)
  })
  deviations.sum <- parLapply(cl, 1:element_num, function(i){
    mod.tmp <- mod78.training.sum[[i]]$mod.tmp
    mu_pred <- predict(mod.tmp, newdata = SCdata.requireZ, what = "mu", type = "response")
    sigma_pred <- predict(mod.tmp, newdata = SCdata.requireZ, what = "sigma", type = "response")
    nu_pred <- predict(mod.tmp, newdata = SCdata.requireZ, what = "nu", type = "response")
    
    dependentvar <- SC_labels[i]
    deviation.df <- data.frame(sub_id=SCdata.requireZ$sub_id)
    observation <- SCdata.requireZ[[dependentvar]]
    centile <- pGG(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred)
    deviation.df[[paste0(SC_labels[i], "_centile")]] <- centile
    deviation.df[[paste0(SC_labels[i], "_deviationZ")]] <- qnorm(centile)
    
    print("breakponit1")
    # compute the standard values while controlling for mean_fd and siteID
    standardvalue.mat <- matrix(NA, nrow(SCdata.requireZ_fixed.F), length(sitelist))
    for (j in 1:length(sitelist)){
      site = sitelist[j]
      SCdata.requireZ_fixed.F$site <- SCdata.requireZ_fixed.M$site <- site
      
      # Female
      mu_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.F, what = "mu", type = "response")
      sigma_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.F, what = "sigma", type = "response")
      nu_pred.fixed.F <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.F, what = "nu", type = "response")
      standardvalue.F <- qGG(centile, mu = mu_pred.fixed.F, sigma = sigma_pred.fixed.F, nu = nu_pred.fixed.F)
      
      # Male
      mu_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.M, what = "mu", type = "response")
      sigma_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.M, what = "sigma", type = "response")
      nu_pred.fixed.M <- predict(mod.tmp, newdata = SCdata.requireZ_fixed.M, what = "nu", type = "response")
      standardvalue.M <- qGG(centile, mu = mu_pred.fixed.M, sigma = sigma_pred.fixed.M, nu = nu_pred.fixed.M)
      
      standardvalue <- (standardvalue.F + standardvalue.M) /2
      standardvalue.mat[,j] <- standardvalue
    }
    print("breakponit2")
    standardvalue <- rowMeans(standardvalue.mat)
    
    deviation.df <- deviation.df %>% drop_na()
    deviation.df[[paste0(SC_labels[i], "_standard")]] <- standardvalue
    
    return(deviation.df)
  })
  saveRDS(deviations.sum, paste0(interfileFolder, "/SCdata.Schaefer200.testset_ADHD_deviation.rds"))
  deviation.requireZ.df <- deviations.sum
}else{
  deviation.requireZ.df <- readRDS(paste0(interfileFolder, "/SCdata.Schaefer200.testset_ADHD_deviation.rds"))

}

deviation.requireZ.df2 <- Reduce(function(x, y) merge(x, y, by = "sub_id", all = TRUE), deviation.requireZ.df)
deviation.requireZ.df3 <- merge(deviation.requireZ.df2, SCdata.test, by="sub_id")
saveRDS(deviation.requireZ.df3, file.path(interfileFolder, "/dataframe_Schaefer200.test.deviation.rds"))

all_colnames <- names(deviation.requireZ.df3)

centile_cols <- all_colnames[grep("_centile$", all_colnames)]
z_cols <- all_colnames[grep("_deviationZ$", all_colnames)]
standard_cols <- all_colnames[grep("_standard$", all_colnames)]
base_cols <- c("sub_id", "age", "sex", "diagnosis", "mean_fd", "TBV", "network_strength_schaefer200", "handedness", "site")


deviation.centile.df <- deviation.requireZ.df3 %>%
  select(all_of(base_cols), all_of(centile_cols))
deviation.centile.df <- deviation.centile.df %>%
  rename_with(~str_remove(.x, "_centile$"), .cols = all_of(centile_cols))
write.csv(deviation.centile.df, 
          file.path(interfileFolder, "dataframe_Schaefer200.test.deviation_centile.csv"), 
          row.names = FALSE)

deviation.Z.df <- deviation.requireZ.df3 %>%
  select(all_of(base_cols), all_of(z_cols))
deviation.Z.df <- deviation.Z.df %>%
  rename_with(~str_remove(.x, "_deviationZ$"), .cols = all_of(z_cols))
write.csv(deviation.Z.df, 
          file.path(interfileFolder, "dataframe_Schaefer200.test.deviation_z.csv"), 
          row.names = FALSE)

deviation.standard.df <- deviation.requireZ.df3 %>%
  select(all_of(base_cols), all_of(standard_cols))
deviation.standard.df <- deviation.standard.df %>%
  rename_with(~str_remove(.x, "_standard$"), .cols = all_of(standard_cols))
write.csv(deviation.standard.df, 
          file.path(interfileFolder, "dataframe_Schaefer200.test.deviation_standard.csv"), 
          row.names = FALSE)


## ADHD
Interest.vars <- c(paste0(SC_labels, "_deviationZ"), "age", "mean_fd")
tableone.df <- CreateTableOne(vars=Interest.vars, strata = "sex", data = deviation.requireZ.df3, test = TRUE)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(interfileFolder, "/demoinfo_sexADHD.csv"), row.names = T)



