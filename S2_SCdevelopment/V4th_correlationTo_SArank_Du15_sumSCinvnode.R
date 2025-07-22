## Validation: large-scale matrix of Du15
#### This script is to conduct correlation analysis 
#### between gam statistical indexes to S-A connectional axis rank.
#### And draw scatter plots & matrix graphs.

library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape) 
library(R.matlab) 


atlas_name <- "Yeo17"
resolution <- 15 
edgenum <- resolution * (resolution + 1) / 2 


homepath <- homepath <- "D:/code/SC_ADHD"
inputFolder <- file.path(homepath, "datasets", "results", "S2") 
functionFolder.SCDev <- file.path(homepath, "code", "functionFolder.SCDev", "gamfunction")
functionFolder.custom <- file.path(homepath, "code", "functions")
FigureFolder <- file.path(homepath, "Figures", atlas_name, "SArank")


if (!dir.exists(FigureFolder)) dir.create(FigureFolder, recursive = TRUE)

source(file.path(functionFolder.SCDev, 'SCrankcorr.R'))
source(file.path(functionFolder.custom, 'plotmatrix.R'))

gamresult.TD.raw <- readRDS(file.path(inputFolder, paste0("gam_stats_TD_res", resolution, ".rds")))
gamresult.ADHD.raw <- readRDS(file.path(inputFolder, paste0("gam_stats_ADHD_res", resolution, ".rds")))


cat("--- 1. Calculating correlation to SC rank ---\n")

clean_gam_results <- function(df) {
    numeric_cols <- c("gam.smooth.F", "gam.smooth.pvalue", "partialRsq", "anova.smooth.zvalue", "anova.smooth.pvalue", "correstimate", "corrp", "gam.edf", "meanderv2")
    cols_to_convert <- intersect(numeric_cols, names(df))
    df[cols_to_convert] <- lapply(df[cols_to_convert], function(x) as.numeric(as.character(x)))
    
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


calculate_sc_rank_correlations <- function(gamresult_df, group_name, resolution) {
    cat(paste("  Processing group:", group_name, "\n"))

    correlation_results <- list()
    
    vars_to_compute <- c("partialRsq2", "meanderv2", "peak.change")
    
    for (computevar in vars_to_compute) {
        if (computevar %in% names(gamresult_df)) {
            correlation_results[[computevar]] <- SCrankcorr(gamresult_df, computevar, resolution, dsdata = FALSE)
        } else {
            cat(paste("    - Variable", computevar, "not found. Skipping.\n"))
        }
    }
    
    final_corr_df <- dplyr::bind_rows(correlation_results)
    print(final_corr_df)
    return(final_corr_df)
}


corr_TD <- calculate_sc_rank_correlations(gamresult.TD, "TD", resolution)
corr_ADHD <- calculate_sc_rank_correlations(gamresult.ADHD, "ADHD", resolution)

# generate scatter plot data
cat("\n--- 2. Generating scatter plots ---\n")

plot_scatter_corr <- function(gamresult_df, dataname, computevar, resolution, FigureFolder) {
    if (!computevar %in% names(gamresult_df)) {
        cat(paste("Cannot create plot for", computevar, "- column not found.\n"))
        return(NULL)
    }
    
    correlation.df <- SCrankcorr(gamresult_df, computevar, resolution, dsdata = TRUE)
    maxth <- max(abs(correlation.df[[computevar]]), na.rm = TRUE)
    
    p <- ggplot(data = correlation.df, aes(x = SCrank, y = .data[[computevar]])) +
        geom_point(aes(color = .data[[computevar]]), size = 5) +
        geom_smooth(method = "lm", color = "black", linewidth = 1.2) +
        scale_color_distiller(type = "seq", palette = "RdBu", direction = -1, limits = c(-maxth, maxth)) +
        scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120)) +
        labs(x = "S-A connectional axis rank", y = computevar) + 
        theme_classic() +
        theme(axis.text = element_text(size = 23, color = "black"),
              axis.title = element_text(size = 23, color = "black"), aspect.ratio = 1,
              axis.line = element_line(linewidth = 0.6),
              axis.ticks = element_line(linewidth = 0.6),
              plot.background = element_rect(fill = "transparent"),
              panel.background = element_rect(fill = "transparent"),
              legend.position = "none")
              
    filename_svg <- file.path(FigureFolder, paste0("scatter_", computevar, "_SCrankcorr_", dataname, ".svg"))
    ggsave(filename_svg, plot = p, dpi = 600, width = 20, height = 14, units = "cm")
    cat(paste("  Scatter plot saved to:", filename_svg, "\n"))
}


for (dataname in c("TD", "ADHD")) {
    gamresult <- if(dataname == "TD") gamresult.TD else gamresult.ADHD
    for (computevar in c("partialRsq2", "meanderv2", "peak.change")) {
        plot_scatter_corr(gamresult, dataname, computevar, resolution, FigureFolder)
    }
}


# plot matrix graphs
cat("\n--- 3. Generating matrix graphs ---\n")

ds.resolution <- resolution
# labels
axeslabels_Yeo17 <- c("VS.P", "SM.A", "VS.C", "SM.B", "DA.A", "DA.B", "FP.C", "DM.c", "VA.A", "TP", "FP.A", "VA.B", "DM.A", "FP.B", "DM.B")

# TD group
dataname_str_TD <- "gamresult.TD"
variable <- "partialRsq2"
lmthr= max(abs(gamresult.TD$partialRsq2), na.rm=T)
Fig.TD.Rsq <- plotmatrix(dataname_str_TD, variable, ds.resolution, Pvar=NA, NAcol="#67001F", lmthr=lmthr, axeslabels_Yeo17, axeslabelsGap=T)
ggsave(file.path(FigureFolder, paste0(atlas_name, "_TD_partialRsq2_matrix.svg")), plot = Fig.TD.Rsq, dpi = 600, width = 20, height = 18, units = "cm")

variable <- "meanderv2"
Fig.TD.derv <- plotmatrix(dataname_str_TD, variable, ds.resolution, Pvar=NA, NAcol="#67001F", lmthr=NA, axeslabels_Yeo17, axeslabelsGap=T)
ggsave(file.path(FigureFolder, paste0(atlas_name, "_TD_meanderv2_matrix.svg")), plot = Fig.TD.derv, dpi = 600, width = 20, height = 18, units = "cm")

# ADHD group
dataname_str_ADHD <- "gamresult.ADHD"
variable <- "partialRsq2"
lmthr= max(abs(gamresult.ADHD$partialRsq2), na.rm=T)
Fig.ADHD.Rsq <- plotmatrix(dataname_str_ADHD, variable, ds.resolution, Pvar=NA, NAcol="#67001F", lmthr=lmthr, axeslabels_Yeo17, axeslabelsGap=T)
ggsave(file.path(FigureFolder, paste0(atlas_name, "_ADHD_partialRsq2_matrix.svg")), plot = Fig.ADHD.Rsq, dpi = 600, width = 20, height = 18, units = "cm")

variable <- "meanderv2"
Fig.ADHD.derv <- plotmatrix(dataname_str_ADHD, variable, ds.resolution, Pvar=NA, NAcol="#67001F", lmthr=NA, axeslabels_Yeo17, axeslabelsGap=T)
ggsave(file.path(FigureFolder, paste0(atlas_name, "_ADHD_meanderv2_matrix.svg")), plot = Fig.ADHD.derv, dpi = 600, width = 20, height = 18, units = "cm")

cat("Matrix plot section is ready. Please uncomment and adapt 'plotmatrix' calls and 'axeslabels' for your Du15 atlas.\n")

print("Analysis and plotting complete.")