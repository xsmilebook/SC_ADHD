library(R.matlab);
library(psych)
library(mgcv);
library(tidyverse)
library(lme4)
library(gamm4)
library(BiocManager)
library(Biostrings)
library(sva)

rm(list = ls())

homepath <- "D:/code/SC_ADHD"

demopath <- paste0(homepath, '/datasets/demography')
interfileFolder <- paste0(homepath, '/datasets/interfileFolder')

functionFolder <- paste0(homepath, "/code/functions")
resultFolder <- paste0(homepath, "/datasets/results")
FigureFolder <- paste0(homepath, "/Figures/Yeo17")

source(paste0(functionFolder, "/ComBat_sva.R"))
Behavior <- read.csv(paste0(demopath, '/THREE_SITES_basic_demo.csv'))



# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
edgenum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2
# load data
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_Yeo', Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))


SCdata$sex <- as.factor(SCdata$sex)
SCdata$diagnosis <- factor(SCdata$ADHD, levels=c("0", "1"), labels=c("TD", "ADHD"))


# Combat
SC_vars <- grep("SC.", names(SCdata), value=T)
model_terms <- c(SC_vars,"scanID", "age", "site","sex", "mean_fd", "diagnosis")
# CV75
comtable <- SCdata %>% select(model_terms) %>%
  drop_na() 
sitetab <- table(comtable$site) #3 sites

batch <- as.factor(comtable$site)
harmonized_data <- data.frame(matrix(NA, nrow(comtable), edgenum))
names(harmonized_data) <- paste0("SC.", c(1:edgenum), "_h")
harmonized_data$scanID <- comtable$scanID

for (i in 1:edgenum){
  ctab <- t(data.matrix(comtable%>%select(SC_vars[i])))
  smooth_var<-"age" ; knots=3; set_fx=TRUE ; covariates="sex+mean_fd+diagnosis"
  region <- SC_vars[i]
  # age model
  modelformula1 <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula1, method="REML", data = comtable)
  mod<- model.matrix(gam.model)
  print(summary(mod))
  # combat. adjust site effect
  combatdata <- ComBat_sva(dat=ctab, batch = batch, mod=mod, par.prior = FALSE)
  harmonized_data[,i]<-t(combatdata)
}

# ## Adopting TD model to remove batch effects of ADHD sample.
# comtable.TD <- comtable %>% filter(diagnosis=="TD")
# for (i in 1:edgenum){
#   ctab <- t(data.matrix(comtable%>%select(SC_vars[i])))
#   smooth_var<-"age" ; knots=3; set_fx=TRUE ; covariates="sex+FD"
#   region <- SC_vars[i]
#   # age model
#   modelformula1 <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
#   gam.model <- gam(modelformula1, method="REML", data = comtable.TD)
#   #mod<- model.matrix(gamm.model$gam, data=comtable)
#   mod <- predict(gam.model, newdata = comtable, type = "lpmatrix")
#   # combat. adjust site effect
#   combatdata <- ComBat_sva(dat=ctab, batch = batch, mod=mod, par.prior = F)
#   harmonized_data[,i]<-t(combatdata)
# }

harmonized_data[,1:edgenum] <- lapply(harmonized_data[,1:edgenum], as.numeric)
dataTable<- merge(harmonized_data, Behavior, by="scanID")
describe(comtable$SC.118)
describe(dataTable$SC.118_h)
corr.test(comtable$SC.118, dataTable$SC.118_h)
saveRDS(dataTable, paste0(interfileFolder, "/combat/SCdata_Yeo", Yeoresolution, "_CV75_sumSCinvnode.sum.msmtcsd.combat_TD_ADHDall_covDiagnose.rds"))
#########################################








