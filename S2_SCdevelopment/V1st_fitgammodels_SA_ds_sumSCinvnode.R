## This script is to fit gam models for each connection for both TD and ADHD groups.

library(mgcv)
library(parallel)
library(tidyverse)


homepath <- "D:/code/SC_ADHD"
interfileFolder <- file.path(homepath, "datasets", "interfileFolder")
resultFolder <- file.path(homepath, "datasets", "results", "S2")
functionFolder.SCDev <- file.path(homepath, "code", "functionFolder.SCDev", "gamfunction")

if (!dir.exists(resultFolder)) {
  dir.create(resultFolder, recursive = TRUE)
}

# Source custom functions
source(file.path(functionFolder.SCDev, "gamsmooth.R"))
source(file.path(functionFolder.SCDev, "plotdata_generate.R"))
cat("Custom functions loaded from:", functionFolder.SCDev, "\n")

# --- 2. Main Analysis Function ---
# This function encapsulates the entire GAM analysis pipeline for a given group.

run_gam_analysis_for_group <- function(group_data, group_name, sc_labels_h, edgenum, smooth_var, covariates, resolution, resultFolder, cl) {
  
  cat(paste0("\n\n========================================================\n"))
  cat(paste0("===== Starting Full Analysis Pipeline for: ", group_name, " =====\n"))
  cat(paste0("========================================================\n\n"))

  # The gam.fit.smooth function uses get(dataname). For mclapply to find the
  # object in forked processes, it must exist in an accessible environment. We
  # temporarily assign the data to a uniquely named object in the global environment.
  dataname <- paste0("current_analysis_data_", group_name)
  assign(dataname, group_data, envir = .GlobalEnv)

  # Also prepare the name for the scaled data object that will be created later.
  dataname2 <- paste0("current_scaled_data_", group_name)
  
  # Use on.exit() to ensure the global variables are removed when the function
  # exits, even if an error occurs.
  on.exit(rm(list = c(dataname, dataname2), envir = .GlobalEnv, inherits = FALSE), add = TRUE)

  # --- Step A: Calculate GAM statistics ---
  cat(paste("Step A: Calculating GAM statistics for", group_name, "group...\n"))
  clusterExport(cl, varlist = c(dataname, "sc_labels_h", "gam.fit.smooth", "smooth_var", "covariates"), envir = environment())

  gam_stats_list <- parLapply(cl, 1:edgenum, function(x) {

    region <- sc_labels_h[x]
    
    gamresult <- gam.fit.smooth(
      region = region,
      dataname = dataname,
      smooth_var = smooth_var,
      covariates = covariates,
      knots = 3,
      set_fx = TRUE,
      stats_only = TRUE,
      mod_only = FALSE
    )
    
    gamresult <- as.data.frame(gamresult)
    gamresult$edge <- region
    return(gamresult)
  })

  # FDR correction
  gam_stats_df <- dplyr::bind_rows(gam_stats_list) # Use bind_rows for safety
  # Ensure columns are numeric before FDR correction
  numeric_cols <- names(gam_stats_df)[sapply(gam_stats_df, is.numeric)]
  gam_stats_df[numeric_cols] <- lapply(gam_stats_df[numeric_cols], as.numeric)
  
  gam_stats_df$pfdr <- p.adjust(gam_stats_df$bootstrap_pvalue, method = "fdr")
  gam_stats_df$sig <- (gam_stats_df$pfdr < 0.05)
  
  stats_filename <- file.path(resultFolder, paste0("gam_stats_", group_name, "_res", resolution, ".rds"))
  saveRDS(gam_stats_df, file = stats_filename)
  cat(paste("Step A complete. GAM statistical results for", group_name, "saved to:", stats_filename, "\n"))

  # --- Step B: Calculate and save full GAM model objects ---
  cat(paste("\nStep B: Calculating and saving full GAM model objects for", group_name, "group...\n"))
  

  gam_models_list <- parLapply(cl, 1:edgenum, function(x) {

    if (x %% 10 == 0) cat("  Fitting model for edge", x, "of", edgenum, "for", group_name, "\n")
    region <- sc_labels_h[x]
    
    model_object <- gam.fit.smooth(
      region = region,
      dataname = dataname,
      smooth_var = smooth_var,
      covariates = covariates,
      knots = 3,
      set_fx = TRUE,
      stats_only = FALSE,
      mod_only = TRUE
    )
    return(model_object)
  })
  
  names(gam_models_list) <- sc_labels_h
  models_filename <- file.path(resultFolder, paste0("gam_models_", group_name, "_res", resolution, "_models.rds"))
  saveRDS(gam_models_list, file = models_filename)
  cat(paste("Step B complete. GAM model objects for", group_name, "saved to:", models_filename, "\n"))

  # --- Step C: Generate plot data and perform scaled analysis ---
  cat(paste("\nStep C: Generating plot data and performing scaled analysis for", group_name, "group...\n"))
  
  agevector <- seq(min(group_data$age, na.rm = TRUE), max(group_data$age, na.rm = TRUE), length.out = 100)
  
  clusterExport(cl, varlist = c("gam_models_list", "plotdata_generate", "agevector"), envir = environment())
  
  plotdatasum <- parLapply(cl, 1:edgenum, function(x) {

    modobj <- gam_models_list[[x]]
    edge_name <- names(gam_models_list)[x]
    
    if (is.null(modobj) || inherits(modobj, "try-error") || all(is.na(modobj))) {
      return(NULL)
    } else {
      plotdata <- plotdata_generate(modobj, dataname, "age", agevector)
      plotdata$SC_label <- edge_name
      return(plotdata)
    }
  })
  
  # KEY FIX: Use dplyr::bind_rows() instead of do.call(rbind, ...).
  # bind_rows is more flexible and handles cases where data frames have
  # different columns by filling missing values with NA.
  plotdatasum.df <- dplyr::bind_rows(plotdatasum)
  
  plot_data_filename <- file.path(resultFolder, paste0("plot_data_", group_name, "_res", resolution, ".rds"))
  saveRDS(plotdatasum.df, file = plot_data_filename)
  cat(paste("Plot data for", group_name, "generated and saved to:", plot_data_filename, "\n"))

  # Scaling analysis
  SCdata.diw <- group_data
  for (i in 1:edgenum) {
    SClabel <- sc_labels_h[i]
    plotdata.tmp <- plotdatasum.df %>% filter(SC_label == SClabel)
    
    if (nrow(plotdata.tmp) > 0) {
      scaling_factor <- plotdata.tmp$fit[1]
      if (!is.na(scaling_factor) && scaling_factor != 0) {
        SCdata.diw[, SClabel] <- SCdata.diw[, SClabel] / scaling_factor
      }
    }
  }
  scaled_data_filename <- file.path(resultFolder, paste0(group_name, "_data_scaled_res", resolution, ".rds"))
  saveRDS(SCdata.diw, file = scaled_data_filename)
  cat(paste("Scaled data for", group_name, "created and saved to:", scaled_data_filename, "\n"))

  # Re-run GAM on scaled data
  cat(paste("Re-running GAM on scaled data for", group_name, "to get final models...\n"))
  
  # As before, assign the scaled data to the globally-accessible object for get() to find
  assign(dataname2, SCdata.diw, envir = .GlobalEnv)
  clusterExport(cl, varlist = c(dataname2), envir = .GlobalEnv)
  
  gam_models_scaled_list <- parLapply(cl, 1:edgenum, function(x) {

    if (x %% 10 == 0) cat("  Fitting scaled model for edge", x, "of", edgenum, "for", group_name, "\n")
    region <- sc_labels_h[x]
    
    model_object <- gam.fit.smooth(
      region = region,
      dataname = dataname2,
      smooth_var = smooth_var,
      covariates = covariates,
      knots = 3,
      set_fx = TRUE,
      stats_only = FALSE,
      mod_only = TRUE
    )
    return(model_object)
  })
  
  names(gam_models_scaled_list) <- sc_labels_h
  scaled_models_filename <- file.path(resultFolder, paste0("gam_models_", group_name, "_scaled_res", resolution, ".rds"))
  saveRDS(gam_models_scaled_list, file = scaled_models_filename)
  
  cat(paste("Step C complete. Final models on scaled data for", group_name, "saved to:", scaled_models_filename, "\n"))
  cat(paste0("===== Finished Analysis for: ", group_name, " =====\n"))
}



resolution <- 15
edgenum <- resolution * (resolution + 1) / 2 # 120
n_cores <- detectCores() - 8
cat(paste("Creating a cluster with", n_cores, "cores...\n"))
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {
  library(tidyr)
  library(mgcv)
  library(gratia)
  library(tidyverse)
  library(dplyr)
  library(ecostats)
})
covariates <- "sex + mean_fd"
smooth_var <- "age"

combat_file <- file.path(interfileFolder, "SCdata_Yeo17_CV75_sumSCinvnode.sum.msmtcsd.combat_match_TD_ADHDall_covDiagnose.rds")
all_data <- readRDS(combat_file)

sc_labels_h <- grep("SC\\..*_h$", names(all_data), value = TRUE)

required_vars <- c("scanID", "age", "sex", "mean_fd", "ADHD", sc_labels_h)
analysis_data <- all_data %>%
  drop_na(all_of(required_vars))

analysis_data$sex <- as.factor(analysis_data$sex)

# Split data into TD and ADHD groups
TD_data <- analysis_data %>% filter(ADHD == 0)
ADHD_data <- analysis_data %>% filter(ADHD == 1)

# TD
run_gam_analysis_for_group(
  group_data = TD_data,
  group_name = "TD",
  sc_labels_h = sc_labels_h,
  edgenum = edgenum,
  smooth_var = smooth_var,
  covariates = covariates,
  resolution = resolution,
  resultFolder = resultFolder,
  cl = cl
)

# ADHD
run_gam_analysis_for_group(
  group_data = ADHD_data,
  group_name = "ADHD",
  sc_labels_h = sc_labels_h,
  edgenum = edgenum,
  smooth_var = smooth_var,
  covariates = covariates,
  resolution = resolution,
  resultFolder = resultFolder,
  cl = cl
)

