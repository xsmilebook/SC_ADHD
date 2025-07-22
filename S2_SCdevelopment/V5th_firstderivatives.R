## This script is to visualize the first derivatives (rate of change)
## of SC strength for both TD and ADHD groups.

library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape)
library(R.matlab)

rm(list = ls())


atlas_name <- "Yeo17"
resolution <- 15
edgenum <- resolution * (resolution + 1) / 2

homepath <- homepath <- "D:/code/SC_ADHD"
stats_inputFolder <- file.path(homepath, "datasets", "results", "S2") 
deriv_inputFolder <- file.path(homepath, "datasets", "results", "S2")

functionFolder.SCDev <- file.path(homepath, "code", "functionFolder.SCDev", "gamfunction")

FigureFolder <- file.path(homepath, "Figures", atlas_name, "deviation")

if (!dir.exists(FigureFolder)) dir.create(FigureFolder, recursive = TRUE)


source(file.path(functionFolder.SCDev, "/SCrankcorr.R")) 

cat("Loading derivative and statistical data...\n")

gamresult.TD.raw <- readRDS(file.path(stats_inputFolder, paste0("gam_stats_TD_res", resolution, ".rds")))
gamresult.ADHD.raw <- readRDS(file.path(stats_inputFolder, paste0("gam_stats_ADHD_res", resolution, ".rds")))

derivativeTD.raw <- readRDS(file.path(deriv_inputFolder, paste0('derivative_df_TD_res', resolution, '.rds')))
derivativeADHD.raw <- readRDS(file.path(deriv_inputFolder, paste0('derivative_df_ADHD_res', resolution, '.rds')))



clean_gam_results <- function(df) {
    numeric_cols <- c("partialRsq") 
    cols_to_convert <- intersect(numeric_cols, names(df))
    if(length(cols_to_convert) > 0) {
        df[cols_to_convert] <- lapply(df[cols_to_convert], function(x) as.numeric(as.character(x)))
    }
    

    if ("partialRsq" %in% names(df)) {
        df <- df %>%
            mutate(
                partialRsq2 = {
                    prsq <- partialRsq
                    mean_prsq <- mean(prsq, na.rm = TRUE)
                    sd_prsq <- sd(prsq, na.rm = TRUE)
                    upper_bound <- mean_prsq + 3 * sd_prsq
                    lower_bound <- mean_prsq - 3 * sd_prsq
                    prsq[prsq > upper_bound | prsq < lower_bound] <- NA
                    prsq
                }
            )
    }
    return(df)
}

gamresult.TD <- clean_gam_results(gamresult.TD.raw)
gamresult.ADHD <- clean_gam_results(gamresult.ADHD.raw)


derivative.dfTD.merge <- derivativeTD.raw %>%
    left_join(dplyr::select(gamresult.TD, edge, partialRsq2), by = c("label_ID" = "edge"))

derivative.dfADHD.merge <- derivativeADHD.raw %>%
    left_join(dplyr::select(gamresult.ADHD, edge, partialRsq2), by = c("label_ID" = "edge"))

cat("Preparing for plotting...\n")


ymin <- min(c(derivative.dfTD.merge$derivative, derivative.dfADHD.merge$derivative), na.rm = TRUE)
ymax <- max(c(derivative.dfTD.merge$derivative, derivative.dfADHD.merge$derivative), na.rm = TRUE)


maxth <- max(abs(c(derivative.dfTD.merge$partialRsq2, derivative.dfADHD.merge$partialRsq2)), na.rm = TRUE)


# TD group
cat("Plotting for TD group...\n")
p_TD <- ggplot(data = derivative.dfTD.merge) +
  geom_line(aes(x = age, y = derivative, group = label_ID, color = partialRsq2), linewidth = 0.8, alpha = 0.9) +
  scale_color_distiller(type = "seq", palette = "RdBu", na.value = "#67001F", limits = c(-maxth, maxth)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) + 
  scale_y_continuous(limits = c(ymin, ymax)) +
  labs(x = "Age (years)", y = "SC Change Rate") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20, color = "black"), 
    axis.title = element_text(size = 20, color = "black"),
    aspect.ratio = 0.9,
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )
  
ggsave(file.path(FigureFolder, paste0(atlas_name, "_FirstDerivative_TD.svg")), plot = p_TD, dpi = 600, width = 13, height = 13, units = "cm")
cat("  TD plot saved.\n")


# ADHD组
cat("Plotting for ADHD group...\n")
p_ADHD <- ggplot(data = derivative.dfADHD.merge) +
  geom_line(aes(x = age, y = derivative, group = label_ID, color = partialRsq2), linewidth = 0.8, alpha = 0.9) +
  scale_color_distiller(type = "seq", palette = "RdBu", na.value = "#67001F", limits = c(-maxth, maxth)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  labs(x = "Age (years)", y = "SC Change Rate") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20, color = "black"), 
    axis.title = element_text(size = 20, color = "black"),
    aspect.ratio = 0.9,
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(FigureFolder, paste0(atlas_name, "_FirstDerivative_ADHD.svg")), plot = p_ADHD, dpi = 600, width = 13, height = 13, units = "cm")
cat("  ADHD plot saved.\n")


print("Analysis and plotting complete.")