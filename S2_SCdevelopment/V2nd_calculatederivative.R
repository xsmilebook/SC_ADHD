## This script generates derivative and posterior derivative values from pre-computed GAM models
## for both TD and ADHD groups in a single run.
Sys.setenv("OPENBLAS_NUM_THREADS" = 1)

library(tidyverse)
library(gratia)
library(mgcv)
library(parallel)
library(MASS) 

rm(list = ls())

resolution <- 15
edgenum <- resolution * (resolution + 1) / 2
n_cores <- detectCores() - 8
cat(paste("Creating a cluster with", n_cores, "cores...\n"))
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {
  library(tidyr)
  library(mgcv)
  library(gratia)
  library(tidyverse)
  library(dplyr)
})


homepath <- "D:/code/SC_ADHD"
inputFolder <- file.path(homepath, "datasets", "results", "S2") 
resultFolder <- file.path(homepath, "datasets", "results", "S2")
functionFolder.SCDev <- file.path(homepath, "code", "functionFolder.SCDev", "gamfunction")

if (!dir.exists(resultFolder)) {
  dir.create(resultFolder, recursive = TRUE)
}

source(file.path(functionFolder.SCDev, 'gamderivatives.R'))
cat("Custom function 'gam.derivatives' loaded.\n")

sc_labels_placeholder <- paste0("SC.", 1:edgenum, "_h")
SCrank <- matrix(NA, resolution, resolution)
for (x in 1:resolution) { for (y in 1:resolution) { SCrank[x,y] = x^2 + y^2 } }
SCrank <- SCrank[lower.tri(SCrank, diag = T)]
SCrank.df <- data.frame(parcel = sc_labels_placeholder, SCrank = SCrank)

groups_to_process <- c("TD", "ADHD")
combat_file <- file.path(homepath, "datasets", "merged_data", "all.SCdata.combat.rds")
all_data <- readRDS(combat_file)
agevector <- all_data$age

clusterExport(cl, varlist = c("gam.derivatives", "agevector", "SCrank.df"))

for (group_name in groups_to_process) {
    
    cat(paste0("\n\n###########################################################\n"))
    cat(paste0("###   Processing Group: ", group_name, "   ###\n"))
    cat(paste0("###########################################################\n"))
    
    cat("--- Loading data and models for", group_name, "group...\n")
    
    gam_models_file <- file.path(inputFolder, paste0("gam_models_", group_name, "_scaled_res", resolution, ".rds"))
    scaled_data_file <- file.path(inputFolder, paste0(group_name, "_data_scaled_res", resolution, ".rds"))

    if (!file.exists(gam_models_file) || !file.exists(scaled_data_file)) {
      warning(paste("Skipping group", group_name, "because input files were not found."))
      next
    }
    
    gammodelsum <- readRDS(gam_models_file)
    group_data_scaled <- readRDS(scaled_data_file)
    

    
    # ***************************************************************

    sc_labels_h <- names(gammodelsum)
    SCrank.df$parcel <- sc_labels_h

    cat("--- Calculating first derivatives for", group_name, "group...\n")
    clusterExport(cl, varlist = c("gammodelsum", "sc_labels_h"))
    derivative.sum <- parLapply(cl, 1:length(gammodelsum), function(x) {
        modobj <- gammodelsum[[x]]
        SClabel.tmp <- sc_labels_h[x]
        if (is.null(modobj)) return(NULL)
        
        tryCatch({
            derivdata <- gam.derivatives(modobj = modobj, 
                                         smooth_var = "age", 
                                         smoothvector = agevector, 
                                         draws = 1, 
                                         increments = 1000, 
                                         return_posterior_derivatives = FALSE)
            
            derivdata$label_ID <- SClabel.tmp
            derivdata$meanSC <- mean(modobj$model[[SClabel.tmp]], na.rm=T)
            return(derivdata)
        }, error = function(e) {
            cat(paste("\nERROR in derivatives for model", SClabel.tmp, ":", e$message, "\n"))
            return(NULL)
        })
    })
    
    derivative.sum_clean <- derivative.sum[!sapply(derivative.sum, is.null)]
    if (length(derivative.sum_clean) == 0) {
        warning(paste("All derivative calculations failed for group", group_name, ". Check errors above."))
        next
    }

    derivative.df <- dplyr::bind_rows(derivative.sum_clean)
    deriv_filename <- file.path(resultFolder, paste0('derivative_df_', group_name, '_res', resolution, '.rds'))
    saveRDS(derivative.df, file = deriv_filename)
    cat(paste("   -> First derivatives saved to:", deriv_filename, "\n"))

    cat("--- Calculating posterior derivatives for", group_name, "group...\n")
    
    derivative.posterior.sum <- parLapply(cl, 1:length(gammodelsum), function(x) {
        modobj <- gammodelsum[[x]]
        SClabel <- sc_labels_h[x]
        if (is.null(modobj)) return(NULL)

        tryCatch({
            derivdata <- gam.derivatives(modobj = modobj, 
                                         smooth_var = "age", 
                                         smoothvector = agevector, 
                                         draws = 1000, 
                                         increments = 1000, 
                                         return_posterior_derivatives = TRUE)
                                         
            derivdata$SCrank <- SCrank.df$SCrank[SCrank.df$parcel == SClabel]
            derivdata$meanSC <- mean(modobj$model[[SClabel]], na.rm=T)
            derivdata$maxSC <- max(modobj$model[[SClabel]], na.rm=T)
            return(derivdata)
        }, error = function(e) {
            cat(paste("\nERROR in posterior derivatives for model", SClabel, ":", e$message, "\n"))
            return(NULL)
        })
    })
    
    non_null_indices <- !sapply(derivative.posterior.sum, is.null)
    derivative.posterior.list <- derivative.posterior.sum[non_null_indices]
    names(derivative.posterior.list) <- sc_labels_h[non_null_indices]
    
    posterior_filename <- file.path(resultFolder, paste0('derivative_posterior_list_', group_name, '_res', resolution, '.rds'))
    saveRDS(derivative.posterior.list, file = posterior_filename)
    cat(paste("   -> Posterior derivatives saved to:", posterior_filename, "\n"))
    
} 

cat("\n===============================================\n")
cat("All derivative analyses finished successfully.\n")
cat("===============================================\n")