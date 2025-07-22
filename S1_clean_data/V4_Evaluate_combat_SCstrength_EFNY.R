library(gridExtra)
library(psych)
library(mgcv);
library(tidyverse)
library(lme4)
library(gamm4)
library(BiocManager)
library(Biostrings)
library(sva)
library(parallel)

rm(list = ls())


Yeoresolution = 17
homepath <- "D:/code/SC_ADHD"

demopath <- paste0(homepath, '/datasets/demography')
interfileFolder <- paste0(homepath, '/datasets/interfileFolder')

functionFolder <- paste0(homepath, "/code/functions")
resultFolder <- paste0(homepath, "/datasets/results")
FigureFolder <- paste0(homepath, "/Figures/Yeo17/combat")

source(paste0(functionFolder, "/ComBat_sva.R"))

# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
edgenum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2
# load data
SCdata.raw <- readRDS(paste0(interfileFolder, '/SCdata_Yeo', Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
SCdata.raw$diagnosis <- factor(SCdata.raw$ADHD, levels=c("0", "1"), labels=c("TD", "ADHD"))
SCdata.raw$sex <- as.factor(SCdata.raw$sex)

combatadd <- "_covDiagnose"
SCdata.combat <- readRDS(paste0(interfileFolder, "/SCdata_Yeo", Yeoresolution, "_CV75_sumSCinvnode.sum.msmtcsd.combat_TD_ADHDall", combatadd,".rds"))
SCdata.combat$diagnosis <- factor(SCdata.combat$ADHD, levels=c("0", "1"), labels=c("TD", "ADHD"))
SCdata.combat$sex <- as.factor(SCdata.combat$sex)

## PCA with density plot per component
SCdata.pca.raw <- principal(SCdata.raw[,2:121], nfactors = 3)
SCdata.pca.combat <- principal(SCdata.combat[,2:121], nfactors = 3)

scores.raw <- as.data.frame(SCdata.pca.raw$scores)
scores.raw$site <- SCdata.raw$site
scores.raw$diagnosis <- SCdata.raw$diagnosis
scores.raw$type <- "raw"

scores.combat <- as.data.frame(SCdata.pca.combat$scores)
scores.combat$site <- SCdata.combat$site
scores.combat$diagnosis <- SCdata.combat$diagnosis
scores.combat$type <- "combat"

scores.merge <- rbind(scores.raw, scores.combat)
scores.merge$diagnosis <- as.factor(scores.merge$diagnosis)
ggplot(data=scores.merge)+
  geom_point(aes(x = RC1, y= RC2, color=site))+
  facet_wrap(~type+diagnosis, ncol=2)+
  theme_minimal()
ggsave(paste0(FigureFolder, "/Combat_test_PCApoint", combatadd,".tiff"), width=20, height = 18, units = "cm")

ggplot(data=scores.merge)+
  geom_density(aes(x = RC1, color=site))+
  facet_wrap(~type+diagnosis, ncol=2)+
  theme_minimal()
ggsave(paste0(FigureFolder, "/Combat_test_PC1density", combatadd,".tiff"), width=20, height = 18, units = "cm")

ggplot(data=scores.merge)+
  geom_density(aes(x = RC2, color=site))+
  facet_wrap(~type+diagnosis, ncol=2)+
  theme_minimal()
ggsave(paste0(FigureFolder, "/Combat_test_PC2density", combatadd,".tiff"), width=20, height = 18, units = "cm")

## Explained variance
Delta_R2 <- data.frame(age = rep(NA, 1), sex =rep(NA, 1), mean_fd=rep(NA, 1), diagnosis=rep(NA, 1), site=rep(NA, 1))
Delta_R2.raw <-Delta_R2
Delta_R2.combat <-Delta_R2

# raw
cl <- makeCluster(8)
clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
invisible(clusterEvalQ(cl, {
  library(gridExtra)
  library(psych)
  library(mgcv);
  library(tidyverse)
  library(lme4)
  library(gamm4)
  library(BiocManager)
  library(Biostrings)
  library(sva)
}))

Delta_R2.rawsum <- parLapply(cl, 1:120, function(x){
  Delta_R2.raw$region <- paste0("SC.", x)
  
  formula0 <- as.formula(sprintf("%s ~ 1", paste0("SC.", x)))
  formula1 <- as.formula(sprintf("%s ~ 1 + age", paste0("SC.", x)))
  formula2 <- as.formula(sprintf("%s ~ 1 + age + sex", paste0("SC.", x)))
  formula3 <- as.formula(sprintf("%s ~ 1 + age + sex + mean_fd", paste0("SC.", x)))
  formula4 <- as.formula(sprintf("%s ~ 1 + age + sex + mean_fd + diagnosis", paste0("SC.", x)))
  formula5 <- as.formula(sprintf("%s ~ 1 + age + sex + mean_fd + diagnosis + site", paste0("SC.", x)))
  
  gam0 <- gam(formula = formula0, data = SCdata.raw, method = "ML")
  Rsq0 <- summary(gam0)$r.sq
  
  gam1 <- gam(formula = formula1, data = SCdata.raw, method = "ML")
  Delta_R2.raw$age[1] <- summary(gam1)$r.sq - Rsq0
  
  gam2 <- gam(formula = formula2, data = SCdata.raw, method = "ML")
  Delta_R2.raw$sex[1] <- summary(gam2)$r.sq - summary(gam1)$r.sq
  
  gam3 <- gam(formula = formula3, data = SCdata.raw, method = "ML")
  Delta_R2.raw$mean_fd[1] <- summary(gam3)$r.sq - summary(gam2)$r.sq
  
  gam4 <- gam(formula = formula4, data = SCdata.raw, method = "ML")
  Delta_R2.raw$diagnosis[1] <- summary(gam4)$r.sq - summary(gam3)$r.sq
  
  gam5 <- gam(formula = formula5, data = SCdata.raw, method = "ML")
  Delta_R2.raw$site[1] <- summary(gam5)$r.sq - summary(gam4)$r.sq
  
  return(Delta_R2.raw)
})

Delta_R2.raw <- do.call(rbind, Delta_R2.rawsum)
Delta_R2.raw[,1:5] <- lapply(Delta_R2.raw[,1:5], as.numeric)


# combat
clusterExport(cl, varlist = ls(), envir = .GlobalEnv)
Delta_R2.combatsum <- parLapply(cl, 1:120, function(x){
  Delta_R2.combat$region <- paste0("SC.", x)
  
  formula0 <- as.formula(sprintf("%s ~ 1", paste0("SC.", x, "_h")))
  formula1 <- as.formula(sprintf("%s ~ 1 + age", paste0("SC.", x, "_h")))
  formula2 <- as.formula(sprintf("%s ~ 1 + age + sex", paste0("SC.", x, "_h")))
  formula3 <- as.formula(sprintf("%s ~ 1 + age + sex + mean_fd", paste0("SC.", x, "_h")))
  formula4 <- as.formula(sprintf("%s ~ 1 + age + sex + mean_fd + diagnosis", paste0("SC.", x, "_h")))
  formula5 <- as.formula(sprintf("%s ~ 1 + age + sex + mean_fd + diagnosis + site", paste0("SC.", x, "_h")))
  
  gam0 <- gam(formula = formula0, data = SCdata.combat, method = "ML")
  Rsq0 <- summary(gam0)$r.sq
  
  gam1 <- gam(formula = formula1, data = SCdata.combat, method = "ML")
  Delta_R2.combat$age[1] <- summary(gam1)$r.sq - Rsq0
  
  gam2 <- gam(formula = formula2, data = SCdata.combat, method = "ML")
  Delta_R2.combat$sex[1] <- summary(gam2)$r.sq - summary(gam1)$r.sq
  
  gam3 <- gam(formula = formula3, data = SCdata.combat, method = "ML")
  Delta_R2.combat$mean_fd[1] <- summary(gam3)$r.sq - summary(gam2)$r.sq
  
  gam4 <- gam(formula = formula4, data = SCdata.combat, method = "ML")
  Delta_R2.combat$diagnosis[1] <- summary(gam4)$r.sq - summary(gam3)$r.sq
  
  gam5 <- gam(formula = formula5, data = SCdata.combat, method = "ML")
  Delta_R2.combat$site[1] <- summary(gam5)$r.sq - summary(gam4)$r.sq
  
  return(Delta_R2.combat)
})

Delta_R2.combat <- do.call(rbind, Delta_R2.combatsum)
Delta_R2.combat[,1:5] <- lapply(Delta_R2.combat[,1:5], as.numeric)

# plot
Delta_R2.raw$Rsq_total <- rowSums(Delta_R2.raw[,1:5])
Delta_R2.raw_long <- Delta_R2.raw %>% pivot_longer(cols=c("age", "sex", "mean_fd", "diagnosis", "site"), names_to="type")
Delta_R2.raw_long$region <- factor(Delta_R2.raw_long$region, levels=Delta_R2.raw$region[order(Delta_R2.raw$Rsq_total, decreasing = T)])
Delta_R2.raw_long$dataproc <- "raw"

Delta_R2.combat$Rsq_total <- rowSums(Delta_R2.combat[,1:5])
Delta_R2.combat_long <- Delta_R2.combat %>% pivot_longer(cols=c("age", "sex", "mean_fd", "diagnosis", "site"), names_to="type")
Delta_R2.combat_long$region <- factor(Delta_R2.combat_long$region, levels=Delta_R2.combat$region[order(Delta_R2.combat$Rsq_total, decreasing = T)])
Delta_R2.combat_long$dataproc <- "combat"

Delta_R2.merge_long <- rbind(Delta_R2.raw_long, Delta_R2.combat_long)
Delta_R2.merge_long$region <- factor(Delta_R2.merge_long$region, levels=Delta_R2.raw$region[order(Delta_R2.raw$Rsq_total, decreasing = T)])
Delta_R2.merge_long$type <- factor(Delta_R2.merge_long$type, levels=rev(c("age", "sex", "mean_fd", "diagnosis", "site")))

ggplot(data = Delta_R2.merge_long, aes(region, y = value, fill = type)) +
  geom_histogram(binwidth = 0.5, stat = "identity",  color = "black", position = "stack",linewidth=0.5) +
  facet_wrap(~dataproc, ncol = 1, nrow = 2, scales = "free") +
  labs(x = "SC edges", y = "R square") +
  scale_x_discrete(breaks = NULL) +
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.5,
        plot.title = element_text(color = "black", size = 14, hjust = 0.5),
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_rect(fill = "transparent",colour = NA),
        axis.title = element_text(color = "black", size = 14),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 14), 
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_text(color = "black", size = 12))


ggsave(paste0(FigureFolder, "/Combat_test", combatadd,".tiff"), width=20, height = 18, units = "cm")



