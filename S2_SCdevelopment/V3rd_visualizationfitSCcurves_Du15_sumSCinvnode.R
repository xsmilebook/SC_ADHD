## This script is to generate fitted values from gam models for both TD and ADHD groups.
## The predicted SC strength will be generated as the age was sampled across the data range, 
## with a total of 1000 data points. Covariates will be set as median or mode.
## The data will be used to draw developmental trajectories.

library(mgcv)
library(tidyverse)
library(parallel)

library(R.matlab)
library(psych)
library(scales)
library(gratia)
library(RColorBrewer)
library(visreg)



atlas_name <- "Yeo17" 
resolution <- 15 
edgenum <- resolution * (resolution + 1) / 2 


homepath <- "D:/code/SC_ADHD"
inputFolder <- file.path(homepath, "datasets", "results", "S2") 

outputFolder <- file.path(homepath, "datasets", "results", "S2")
functionFolder.SCDev <- file.path(homepath, "code", "functionFolder.SCDev", "gamfunction")
FigureFolder <- file.path(homepath, "Figures", "Yeo17", "fitSCcurves")


if (!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
if (!dir.exists(FigureFolder)) dir.create(FigureFolder, recursive = TRUE)


source(file.path(functionFolder.SCDev, 'plotdata_generate.R'))
# source(file.path(functionFolder.SCDev, '/plotdata_derivatives.R'))
# source(file.path(functionFolder.SCDev, '/gammsmooth.R'))
cat("Custom functions loaded.\n")

cat("Loading pre-computed models and statistical results...\n")

gammodelsumTD <- readRDS(file.path(inputFolder, paste0("gam_models_TD_scaled_res", resolution, ".rds")))
gamresultsumTD <- readRDS(file.path(inputFolder, paste0("gam_stats_TD_res", resolution, ".rds")))

gammodelsumADHD <- readRDS(file.path(inputFolder, paste0("gam_models_ADHD_scaled_res", resolution, ".rds")))
gamresultsumADHD <- readRDS(file.path(inputFolder, paste0("gam_stats_ADHD_res", resolution, ".rds")))


analysis_data_file <- file.path(homepath, "datasets", "merged_data", "all.SCdata.combat.rds")
analysis_data <- readRDS(analysis_data_file)

analysis_data <- analysis_data %>% 
    drop_na(all_of(c("age", "sex", "mean_fd"))) %>%
    mutate(sex = as.factor(sex))

agevector <- analysis_data$age

# --- 6. SCrank
MatrixSCrank <- matrix(NA, nrow=resolution, ncol=resolution)
indexup <- upper.tri(MatrixSCrank)
indexsave <- !indexup
MatrixSCrank.index <- MatrixSCrank
MatrixSCrank.index[indexsave] <- 1:edgenum

for (x in 1:resolution) {
  for (y in 1:resolution) {
    MatrixSCrank[x, y] <- (x + y)^2 + (x - y)^2
  }
}
MatrixSCrank[indexup] <- NA

ranked_values <- rank(MatrixSCrank[indexsave], ties.method = "first")


gamresultsumTD$SCrank <- ranked_values
gamresultsumADHD$SCrank <- ranked_values

# --- 7. fitting data generation
cat("Generating fitted values data for both groups...\n")

groups <- c("TD", "ADHD")
plotdatasum <- list()

for (group_name in groups) {
    
    cat(paste("  Processing group:", group_name, "\n"))
    
    gammodelsum <- if (group_name == "TD") gammodelsumTD else gammodelsumADHD
    gamresultsum <- if (group_name == "TD") gamresultsumTD else gamresultsumADHD
    

    dataname <- paste0("temp_data_for_", group_name)
    group_data <- if(group_name == "TD") filter(analysis_data, ADHD == 0) else filter(analysis_data, ADHD == 1)
    assign(dataname, group_data, envir = .GlobalEnv)
    
    plot_data_list <- list()
    for (x in 1:edgenum) {
        modobj <- gammodelsum[[x]]
        if (is.null(modobj)) next # 如果模型为空则跳过
        
        plotdata <- plotdata_generate(modobj, dataname = dataname, smooth_var = "age", smoothvector = agevector)
        
        plotdata$SC_label <- names(modobj$model)[1] 
        plotdata$SCrank <- gamresultsum$SCrank[x]
        plotdata$PartialRsq <- gamresultsum$partialRsq[x] 
        # plotdata$meanderv2 <- gamresultsum$meanderv2[x]
        
        plot_data_list[[x]] <- plotdata
    }
    
    rm(list = dataname, envir = .GlobalEnv)
    
    plotdatasum.df <- dplyr::bind_rows(plot_data_list)
    plotdatasum[[group_name]] <- plotdatasum.df 
    
    saveRDS(plotdatasum.df, file.path(outputFolder, paste0("plotdatasum_scaled_", atlas_name, "_", group_name, ".rds")))
}

# --- 8. plotting developmental trajectories ---
cat("Plotting developmental trajectories...\n")

plotdatasum.df.TD <- plotdatasum[["TD"]]
plotdatasum.df.ADHD <- plotdatasum[["ADHD"]]

ymin <- min(c(plotdatasum.df.TD$fit, plotdatasum.df.ADHD$fit), na.rm = TRUE)
ymax <- max(c(plotdatasum.df.TD$fit, plotdatasum.df.ADHD$fit), na.rm = TRUE)


gamresultsumTD <- gamresultsumTD %>%
  mutate(

    partialRsq = as.numeric(as.character(partialRsq)),
    partialRsq2 = if ("partialRsq" %in% names(.)) {
        ifelse(partialRsq > mean(partialRsq, na.rm=T) + 3*sd(partialRsq, na.rm=T), NA, partialRsq)
    } else { NA } 
  )
maxth <- max(abs(gamresultsumTD$partialRsq2), na.rm=T)
plotdatasum.df.TD <- plotdatasum.df.TD %>% left_join(dplyr::select(gamresultsumTD, edge, partialRsq2), by = c("SC_label" = "edge"))

p_TD <- ggplot(plotdatasum.df.TD, aes(x=age, y=fit, group=SC_label)) +
  geom_line(aes(color=partialRsq2), linewidth=0.8, alpha=0.8) +
  scale_color_distiller(type="seq", palette = "RdBu", limits=c(-maxth, maxth), na.value="#67001F", direction = -1) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  labs(x="Age (years)", y="SC strength (ratio)", title = "TD Group Developmental Trajectories") +
  theme_classic() +
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20, color="black"), aspect.ratio = 0.9,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=20, hjust = 0.5), legend.position = "none")

ggsave(file.path(FigureFolder, paste0(atlas_name, "_devcurve_TD.svg")), plot = p_TD, dpi=600, width=14, height=13, units="cm")
cat("TD plot saved.\n")

# --- ADHD
gamresultsumADHD <- gamresultsumADHD %>%
  mutate(
    partialRsq = as.numeric(as.character(partialRsq)),
    partialRsq2 = if ("partialRsq" %in% names(.)) {
        ifelse(partialRsq > mean(partialRsq, na.rm=T) + 3*sd(partialRsq, na.rm=T), NA, partialRsq)
    } else { NA }
  )
plotdatasum.df.ADHD <- plotdatasum.df.ADHD %>% left_join(dplyr::select(gamresultsumADHD, edge, partialRsq2), by = c("SC_label" = "edge"))

p_ADHD <- ggplot(plotdatasum.df.ADHD, aes(x=age, y=fit, group=SC_label)) +
  geom_line(aes(color=partialRsq2), linewidth=0.8, alpha=0.8) +
  scale_color_distiller(type="seq", palette = "RdBu", limits=c(-maxth, maxth), na.value="#67001F", direction = -1) +
  scale_y_continuous(limits = c(ymin, ymax)) +
  labs(x="Age (years)", y="SC strength (ratio)", title = "ADHD Group Developmental Trajectories") +
  theme_classic() +
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20, color="black"), aspect.ratio = 0.9,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=20, hjust = 0.5), legend.position = "none")

ggsave(file.path(FigureFolder, paste0(atlas_name, "_devcurve_ADHD.svg")), plot = p_ADHD, dpi=600, width=14, height=13, units="cm")
cat("ADHD plot saved.\n")

print("Analysis and plotting complete.")