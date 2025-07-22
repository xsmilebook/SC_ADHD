library(R.matlab)
library(psych)
library(mgcv)     
library(tidyverse)
# library(lme4)
# library(gamm4)
library(BiocManager)
library(Biostrings)
library(sva)


rm(list = ls())


datapath <- 'D:/code/R_projects/SC_ADHD'
interfileFolder <- paste0(datapath, '/datasets/combine_folder')
homepath <- "D:/code/SC_ADHD"
functionFolder <- paste0(homepath, "/code/functions")
resultFolder <- paste0(datapath, "/datasets/combine_folder/combat")
source(paste0(functionFolder, "/ComBat_sva.R"))

resolution <- 15
edgenum <- resolution * (resolution + 1) / 2 # 15*16/2 = 120

SCdata <- readRDS(paste0(interfileFolder, '/all.SCdata_demo.sum.msmtcsd.rds'))

SCdata$sex <- factor(SCdata$sex, levels = c("M", "F"))


SCdata <- SCdata %>%
  mutate(
    diagnosis_label = factor(ADHD, levels = c(0, 1), labels = c("TD", "ADHD"))
  )

Behavior <- SCdata %>% select(!starts_with("SC."))
if (!"subID" %in% names(Behavior)) {
  stop("Error: 'subID' column not found in the data. Cannot proceed without a unique identifier.")
}


SC_vars <- grep("SC.", names(SCdata), value = TRUE)
SCdata <- SCdata %>% mutate(WB_SCmean = rowMeans(select(., all_of(SC_vars))))
n0 <- nrow(SCdata)
boxplot(SCdata$WB_SCmean, main="Distribution of WB_SCmean before cleaning")

SCdata <- SCdata %>%
  filter(WB_SCmean > mean(WB_SCmean) - 3 * sd(WB_SCmean),
         WB_SCmean < mean(WB_SCmean) + 3 * sd(WB_SCmean))
n1 <- nrow(SCdata)
print(paste(n1, "obs included,", (n0 - n1), "obs were removed due to extreme deviation of WB_SCmean."))



model_terms <- c(SC_vars, "subID", "age", "site", "sex", "mean_fd", "diagnosis_label")

comtable <- SCdata %>% select(all_of(model_terms)) %>% drop_na()

print("Sample size per site:")
print(table(comtable$site))

batch <- as.factor(comtable$site)
harmonized_data <- data.frame(matrix(NA, nrow(comtable), edgenum))
names(harmonized_data) <- paste0("SC.", 1:edgenum, "_h")
harmonized_data$subID <- comtable$subID

for (i in 1:edgenum) {
  if (i %% 10 == 0) {
    print(paste("Processing edge", i, "of", edgenum, "..."))
  }
  
  ctab <- t(data.matrix(comtable %>% select(all_of(SC_vars[i]))))
  
  smooth_var <- "age"
  knots <- 3
  set_fx <- TRUE
  covariates <- "sex+mean_fd+diagnosis_label" 
  region <- SC_vars[i]
  
  
  modelformula1 <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  
  gam.model <- tryCatch({
    gam(modelformula1, data = comtable)
  }, error = function(e) {
    print(paste("Error in GAM for edge", region, ":", e$message))
    return(NULL)
  })
  
  if (is.null(gam.model)) {
    harmonized_data[, i] <- NA
    next
  }
  
  mod <- model.matrix(gam.model)
  
  combatdata <- ComBat_sva(dat = ctab, batch = batch, mod = mod, par.prior = FALSE)
  harmonized_data[, i] <- t(combatdata)
}

harmonized_data[, 1:edgenum] <- lapply(harmonized_data[, 1:edgenum], as.numeric)

dataTable <- merge(harmonized_data, Behavior, by = "subID")

print("--- ComBat Sanity Check ---")
check_edge <- if("SC.21" %in% SC_vars) "SC.21" else SC_vars[1]
check_edge_h <- paste0(check_edge, "_h")

print(paste("Checking edge:", check_edge))
print("Original data description:")
print(describe(comtable[[check_edge]]))
print("Harmonized data description:")
print(describe(dataTable[[check_edge_h]]))
print("Correlation between original and harmonized:")
corr.test(comtable[[check_edge]], dataTable[[check_edge_h]], use = "pairwise.complete.obs")

saveRDS(dataTable, paste0(resultFolder, "/all.SCdata.combat.rds"))

print("Analysis complete. Harmonized data saved.")