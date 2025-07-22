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
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
element_num <- Yeoresolution.delLM*(Yeoresolution.delLM+1)/2
# input directory
homepath <- "D:/code/SC_ADHD"

demopath <- paste0(homepath, '/demography')
interfileFolder <- file.path(homepath, "datasets", "interfileFolder")
functionFolder <- paste0(homepath, "/code/functions")
resultFolder <- paste0(homepath, "/datasets/results/S30")
functionFolder.SCDev <- paste0(homepath, "/code/functionFolder.SCDev/gamfunction")
FigureFolder <- paste0(homepath, '/Figures/Yeo17/constructNM')

# Load data
SCdata <- readRDS(file.path(interfileFolder, "SCdata_Yeo17_CV75_sumSCinvnode.sum.msmtcsd.merge.rds"))
SCdata$sex <- as.factor(SCdata$sex)
SCdata <- SCdata %>% rename(siteID=site)
## Define diagnosis variables
SCdata <- SCdata %>% mutate(
  diagnosis_label=case_when(
    ADHD==1 ~ "ADHD",
    ADHD==0 ~ "TD",
    .default = NA
  )
)

SCdata$diagnosis_label <- factor(SCdata$diagnosis_label, levels = c("TD", "ADHD"))


# data clean
SCdata <- SCdata %>% mutate(WB_SCmean = rowMeans(select(.,starts_with("SC."))))
n0 = nrow(SCdata)
boxplot(SCdata$WB_SCmean)
SCdata <- SCdata %>% filter(WB_SCmean>mean(WB_SCmean)-3*sd(WB_SCmean), WB_SCmean<mean(WB_SCmean)+3*sd(WB_SCmean))
n1 = nrow(SCdata)
print(paste(n1, "obs included,", (n0-n1), "obs were removed due to extreme deviation of WB_SCmean."))


SCdata.ADHDTD <- SCdata %>% filter(diagnosis_label=="ADHD" | diagnosis_label=="TD")
SCdata.TD <- SCdata.ADHDTD %>% filter(diagnosis_label=="TD")
SCdata.ADHD <- SCdata.ADHDTD %>% filter(diagnosis_label=="ADHD")

# source function
source(paste0(functionFolder, "/Construct_gamlss_set.R"))

SCdata.TD.trainset <- readRDS(paste0(resultFolder, "/SCdata.TD.trainset_SCYeo", element_num, ".rds"))
SCdata.TD.testset <- readRDS(paste0(resultFolder, "/SCdata.TD.testset_SCYeo", element_num, ".rds"))
SCdata.TD.trainset <- SCdata.TD.trainset %>% mutate(
  diagnosis_label=case_when(
    ADHD==1 ~ "ADHD",
    ADHD==0 ~ "TD",
    .default = NA
  )
)

SCdata.TD.trainset$diagnosis_label <- factor(SCdata.TD.trainset$diagnosis_label, levels = c("TD", "ADHD"))

SCdata.TD.testset <- SCdata.TD.testset %>% mutate(
  diagnosis_label=case_when(
    ADHD==1 ~ "ADHD",
    ADHD==0 ~ "TD",
    .default = NA
  )
)

SCdata.TD.testset$diagnosis_label <- factor(SCdata.TD.testset$diagnosis_label, levels = c("TD", "ADHD"))

# Describe demographic information.
SCdata.ADHDTD <- SCdata.ADHDTD %>% mutate(
  group = case_when(
    scanID %in% SCdata.TD.trainset$scanID ~ "TDtrain",
    scanID %in% SCdata.TD.testset$scanID ~ "TDtest",
    ADHD == 1 ~ "ADHD",
    .default = NA
  ))

table(SCdata.ADHDTD$group)
demovar <- c("sex", "age", "mean_fd", "siteID")
tableone.df <- CreateTableOne(demovar, strata="group", data=SCdata.ADHDTD, 
                             factorVars = c("sex", "siteID"), includeNA = T,
                             test=T)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_testtrainsets_Yeo17.csv"), row.names = T)

tableone.df <- CreateTableOne(demovar, strata="group", data=SCdata.ADHDTD[SCdata.ADHDTD$diagnosis_label=="TD", ], 
                              factorVars = c("sex", "siteID"), includeNA = T,
                              test=T)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_testtrainsetsTD_Yeo17.csv"), row.names = T)

## Fit normative models
SCdata.TD.trainset <- as.data.frame(SCdata.TD.trainset)
SCdata.TD.trainset[,c("sex", "siteID")] <- lapply(SCdata.TD.trainset[,c("sex", "siteID")], as.factor)
dataname <- "SCdata.TD.trainset"
smoothterm <- "age"
covariates <- "sex+mean_fd"
randomvar <- "siteID"
mu.df <- sigma.df <- degree <- 2
distribution.fam <- "JSU"
IDvar <- "scanID"
quantile.vec <- c(0.025, 0.5, 0.975)
stratify <- c("sex", "siteID")

# MARK: debug


if (! file.exists(paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,".TDtraining.sum.rds"))){
  cl <- makeCluster(30)
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
      dependentvar <- paste0("SC.", i)
      sumlist <- construct_gamlss(dataname, dependentvar, smoothterm, covariates, randomvar, 
                                  mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify)
      
      return(sumlist)
    }, error = function(e){

      error_info <- list(
        status = "error",
        index = i,
        dependent_variable = paste0("SC.", i),
        error_message = e$message # 捕获具体的错误信息
      )
      
      
      return(error_info)
    })
  })
  saveRDS(mod78.training.sum, paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,".TDtraining.sum.rds"))
}else{
  mod78.training.sum <- readRDS(paste0(interfileFolder, "/GAMLSS_Yeo", Yeoresolution,".TDtraining.sum.rds"))
}

# compute deviations
SCdata.TD.trainset.subset <- SCdata.TD.trainset %>% select(c(paste0("SC.", 1:element_num), "age", "sex", "siteID","scanID",
                                                                                     "mean_fd", "diagnosis_label")) %>% drop_na()
SCdata.ADHD <- SCdata.ADHD %>% filter(siteID %in% SCdata.TD.trainset$siteID)
SCdata.ADHD.subset <- SCdata.ADHD %>% select(c(paste0("SC.", 1:element_num), "age", "sex", "siteID","scanID",
                                                                       "mean_fd","diagnosis_label")) %>% drop_na()



SCdata.TD.testset.subset <- SCdata.TD.testset %>% select(c(paste0("SC.", 1:element_num), "age", "sex", "siteID","scanID",
                                                                                   "mean_fd", "diagnosis_label")) %>% drop_na()
SCdata.requireZ <- rbind(SCdata.ADHD.subset, SCdata.TD.testset.subset)

# compute deviation & draw plots
SCdata.requireZ$sex <- as.factor(SCdata.requireZ$sex)
SCdata.requireZ_fixed <- SCdata.requireZ
SCdata.requireZ_fixed$mean_fd <- mean(SCdata.requireZ_fixed$mean_fd)
SCdata.requireZ_fixed.F <- SCdata.requireZ_fixed.M <- SCdata.requireZ_fixed

SCdata.TD.trainset$siteID <- droplevels(SCdata.TD.trainset$siteID)
sitelist <- unique(SCdata.TD.trainset$siteID)

gam.data2 <- SCdata.TD.trainset.subset






if (! file.exists(paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".testset_ADHD_deviation.rds"))){
  cl <- makeCluster(30)
  
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
    
    dependentvar <- paste0("SC.", i)
    deviation.df <- data.frame(scanID=SCdata.requireZ$scanID)
    observation <- SCdata.requireZ[[dependentvar]]
    centile <- pGG(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred)
    deviation.df[[paste0("SC.", i, "_centile")]] <- centile
    deviation.df[[paste0("SC.", i, "_deviationZ")]] <- qnorm(centile)
    
    print("breakponit1")
    # compute the standard values while controlling for mean_fd and siteID
    standardvalue.mat <- matrix(NA, nrow(SCdata.requireZ_fixed.F), length(sitelist))
    for (j in 1:length(sitelist)){
      siteID = sitelist[j]
      SCdata.requireZ_fixed.F$siteID <- SCdata.requireZ_fixed.M$siteID <- siteID
      
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
    deviation.df[[paste0("SC.", i, "_standard")]] <- standardvalue
    
    return(deviation.df)
  })
  saveRDS(deviations.sum, paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".testset_ADHD_deviation.rds"))
  deviation.requireZ.df <- deviations.sum
}else{
  deviation.requireZ.df <- readRDS(paste0(interfileFolder, "/SCdata.sum75_Yeo", Yeoresolution,".testset_ADHD_deviation.rds"))

}
deviation.requireZ.unique <- lapply(deviation.requireZ.df, function(x) x %>% distinct(scanID, .keep_all = T))

deviation.requireZ.df2 <- Reduce(function(x, y) merge(x, y, by = "scanID", all = TRUE), deviation.requireZ.unique)
deviation.requireZ.df3 <- merge(deviation.requireZ.df2, SCdata.ADHDTD, by="scanID") %>% 
  filter(scanID %in% c(SCdata.ADHD$scanID, SCdata.TD.testset$scanID))
deviation.requireZ.df3$ADHD <- factor(deviation.requireZ.df3$ADHD, levels = c(0, 1), labels = c("TD", "ADHD"))
describeBy(deviation.requireZ.df3$SC.100_deviationZ, group=deviation.requireZ.df3$if_TD)

saveRDS(deviation.requireZ.df3, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.deviations_TDtest_ADHD.rds'))



## ADHD
deviation.requireZ.ADHD <- deviation.requireZ.df3 %>% filter(diagnosis_label=="ADHD")
Interest.vars <- c(paste0("SC.", 1:element_num, "_deviationZ"), "age", "mean_fd")
tableone.df <- CreateTableOne(vars=Interest.vars, strata = "sex", data = deviation.requireZ.ADHD, test = TRUE)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_sexADHD.csv"), row.names = T)
## TD
deviation.requireZ.TD <- deviation.requireZ.df3 %>% filter(diagnosis_label=="TD")
Interest.vars <- c(paste0("SC.", 1:element_num, "_deviationZ"), "age", "mean_fd")
tableone.df <- CreateTableOne(vars=Interest.vars, strata = "sex", data = deviation.requireZ.TD, test = TRUE)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_sexTD.csv"), row.names = T)


