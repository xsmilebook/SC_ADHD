library(MatchIt)
library(tableone)
library(dplyr)


homepath <- "D:/code/SC_ADHD"
interfileFolder <- file.path(homepath, "datasets", "interfileFolder")
resultFolder <- file.path(homepath, "datasets", "results", "S2")
combat_file <- file.path(interfileFolder, "SCdata_Yeo17_CV75_sumSCinvnode.sum.msmtcsd.combat_TD_ADHDall_covDiagnose.rds")

SCdata.base <- readRDS(combat_file)

SCdata.base$if_TD <- as.integer(!(SCdata.base$ADHD == 1))

SCdata.base.m <- matchit(if_TD~sex + age + mean_fd, data = SCdata.base, method = "nearest", ratio=1.5, caliper = 0.02)
summary(SCdata.base.m)

SCdata.base.match <- match.data(SCdata.base.m)
# filter the satisfied TD subject
SCdata.base.final <- SCdata.base %>% 
  filter(!(if_TD == 1 & !scanID %in% SCdata.base.match$scanID))



saveRDS(SCdata.base.final, paste0(interfileFolder, "/SCdata_Yeo17_CV75_sumSCinvnode.sum.msmtcsd.combat_match_TD_ADHDall_covDiagnose.rds"))

# description
demovar <- c("sex", "age", "mean_fd", "site")

tableone.df <- CreateTableOne(demovar, strata="if_TD", data=SCdata.base.match, 
                              factorVars = c("sex", "site"), includeNA = T,
                              test=T)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder, "/demoinfo_ADHD_TD_matched_Yeo17.csv"), row.names = T)