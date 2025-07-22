## THREE SITES data（CCNP/EFNY/PKU6）
## This script is to generate a dataframe, in which each column is the strength for an edge.
## For schaefer 400 atlas, 70786 edges left after deleting edges connecting to limbic regions.
library(R.matlab)
library(tidyverse)
library(parallel)
library(openxlsx)
library(rjson)
library(corrplot)
rm(list = ls())
# Set atlas
Yeoresolution <- 17
# Set path and load data

homepath <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD"
SC_path_EFNY <-'/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/defaultatlas'
SC_path_PKU6 <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/PKU6/SCmat'
SC_path_CCNP <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/SC'
Volume_path_EFNY <-'/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/schaefer400_nodevolume'
Volume_path_PKU6 <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/PKU6/schaefer400_nodevolume'
Volume_path_CCNP <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/schaefer400_nodevolume'

demopath <- paste0(homepath, '/datasets/demography')
atlasFolder <- paste0(homepath, '/datasets/atlas/schaefer400')
interfileFolder <- paste0(homepath, '/datasets/interfileFolder')
functionFolder <- paste0(homepath, "/code/functions")
resultFolder <- paste0(homepath, "/datasets/results")
FigureFolder <- paste0(homepath, "/Figures/Yeo17")

Behavior <- read.csv(paste0(demopath, '/THREE_SITES_basic_demo.csv'))

#### import schaefer400 index
schaefer400_index_SA<-read.csv(paste0(atlasFolder, '/schaefer400_index_SA.csv'))

## qsiprep output matrix is in Yeo 7 order, so reorder schaefer400 index to Yeo 7 order
schaefer400_index<-schaefer400_index_SA[order(schaefer400_index_SA$index),]
limbicindex <- which(str_detect(schaefer400_index$label_17network, "Limbic"))
schaefer400_index <- schaefer400_index[-limbicindex, ]
schaefer376_delLM <- schaefer400_index$index
## Rearrange left and right regions.
schaefer400_index$index_7network_LRmixed <- schaefer400_index$index
schaefer400_index$index_7network_LRmixed[str_detect(schaefer400_index$label, "RH_")] <- schaefer400_index$index_7network_LRmixed[str_detect(schaefer400_index$label, "RH_")] - 200
orderYeo_7<-order(schaefer400_index$index_7network_LRmixed)

schaefer400_index$index_17network_LRmixed <- schaefer400_index$index_17network
schaefer400_index$index_17network_LRmixed[str_detect(schaefer400_index$label_17network, "RH_")] <- schaefer400_index$index_17network_LRmixed[str_detect(schaefer400_index$label_17network, "RH_")] - 200
orderYeo_17<-order(schaefer400_index$index_17network_LRmixed)

#### import SC data
#### 376 regions, 377*376/2=70876 SCs
#################################################
colname <- character(length = 70876)
for (i in 1:70876){
  colname[i] <- paste0('SC.', as.character(i))
}

SCdata.sum <- mclapply(1:nrow(Behavior), function(i){
  SCdat <- NULL
  scanID <- Behavior$scanID[i]
  site <- Behavior$site[i]
  
  if (site == "EFNY"){
    SCname <- paste0(scanID, '_dir-PA_space-T1w_desc-preproc_msmtconnectome.mat')
    SC_file_path <- paste0(SC_path_EFNY, '/', SCname)
  }else if (site == "PKU6"){
    SCname <- paste0(scanID, '.mat')
    SC_file_path <- paste0(SC_path_PKU6, '/', SCname)
  }else{
    SCname <- paste0(scanID, '_space-T1w_desc-preproc_dhollanderconnectome.mat')
    SC_file_path <- file.path(SC_path_CCNP, SCname)
  }
  
  if (file.exists(SC_file_path)){
    SCmat <- readMat(SC_file_path)
    SCmat <- SCmat$schaefer400.sift.invnodevol.radius2.count.connectivity[schaefer376_delLM, schaefer376_delLM]
    if (Yeoresolution == 7){
      SCmat <- SCmat[orderYeo_7, orderYeo_7]
    }else if (Yeoresolution == 17){
      SCmat <- SCmat[orderYeo_17, orderYeo_17]
    }else{
      print("Invalid Yeoresolution!")
    }
    
    indexup <- upper.tri(SCmat)
    indexsave <- !indexup ###keep lower triangle and diagonal
    SCdat <- as.data.frame(c(SCmat[indexsave]))
    SCdat <- as.data.frame(t(SCdat), row.names = NULL)
    names(SCdat) <- colname
    row.names(SCdat) <- NULL
    SCdat$scanID[1] <- scanID
  }else {
    print(paste("File not found for ID:", scanID, "at path:", SC_file_path))
  }
  return(SCdat)
}, mc.cores = 50)
ncoldf <- lapply(SCdata.sum, function(x) ncol(x))
SCdata.df <- do.call(rbind, SCdata.sum)
saveRDS(SCdata.df, paste0(interfileFolder, '/SCdataYeo', Yeoresolution,'.sum.msmtcsd.delLM.rds'))
SCdata.sum.merge <- merge(SCdata.df, Behavior, by="scanID")
## calculate CV
meanSC<-colMeans(SCdata.sum.merge[,2:70877])
sd.SC <- mclapply(1:70876, function(x) {
  sd.tmp<-sd(SCdata.sum.merge[,x+1])
  return(sd.tmp)
}, mc.cores = 40)
sd.SC<-as.numeric(sd.SC)
CV.SC<-sd.SC/meanSC
Perct.CV.SC <- quantile(CV.SC, probs=seq(0, 1, 0.25)) # extract 
# 0%         25%        50%        75%       100% 
# 0.3510151  0.9908105  1.2180421  1.5238890 11.0477614
## ID of edges over threshold
deleteindex.delLM <- which(CV.SC>Perct.CV.SC[4])
SCdata.sum.merge[,deleteindex.delLM+1] <- 0
meanSC[deleteindex.delLM] <-0
saveRDS(SCdata.sum.merge, paste0(interfileFolder, '/SCdata.sum.CV75.merge.Yeo', Yeoresolution,'.delLM.rds'))
saveRDS(deleteindex.delLM, paste0(interfileFolder, '/CV75_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))

# Validation Threshold = 25th CV
deleteindex.delLM.25 <- which(CV.SC>Perct.CV.SC[2])
SCdata.sum.merge.CV25 <- SCdata.sum.merge
SCdata.sum.merge.CV25[,deleteindex.delLM.25+1] <-0
meanSC.CV25 <- meanSC; meanSC.CV25[deleteindex.delLM.25] <-0
saveRDS(SCdata.sum.merge.CV25, paste0(interfileFolder, '/SCdata.sum.CV25.merge.Yeo', Yeoresolution,'.delLM.rds'))
saveRDS(deleteindex.delLM.25, paste0(interfileFolder, '/CV25_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))
########################################################################

## plot
SCdata.sum.merge <- readRDS(paste0(interfileFolder, '/SCdata.sum.CV75.merge.Yeo', Yeoresolution,'.delLM.rds'))
SCdata.sum.merge.CV25 <- readRDS(paste0(interfileFolder, '/SCdata.sum.CV25.merge.Yeo', Yeoresolution,'.delLM.rds'))
meanSC<-colMeans(SCdata.sum.merge[,2:70877])
meanSC25<-colMeans(SCdata.sum.merge.CV25[,2:70877])
#376
Matsize<-376
Matrix.376 <- matrix(NA, nrow=Matsize, ncol =Matsize)
indexup <- upper.tri(Matrix.376)
indexsave <- !indexup ###keep lower triangle and diagonal
index <- as.numeric(meanSC)
Matrix.376[indexsave] <- index
Matrix.376[indexup] <- t(Matrix.376)[indexup]
colnames(Matrix.376) <-seq(1, Matsize)
rownames(Matrix.376) <-seq(1, Matsize)
tiff( 
  filename = paste0(FigureFolder, '/SCmatrix/SClog376_CV75_Yeo', Yeoresolution,'.tiff'),
  width = 600, 
  height = 600,
  units = "px",
  bg = "white",
  res = 100)
image(log(Matrix.376), col=rev(COL2(diverging = "RdBu", n=200)), axes = FALSE)
dev.off()

# CV25
index <- as.numeric(meanSC25)
Matrix.376[indexsave] <- index
Matrix.376[indexup] <- t(Matrix.376)[indexup]
colnames(Matrix.376) <-seq(1, Matsize)
rownames(Matrix.376) <-seq(1, Matsize)
tiff( 
  filename = paste0(FigureFolder, '/SCmatrix/SClog376_CV25_Yeo', Yeoresolution,'.tiff'),
  width = 600, 
  height = 600,
  units = "px",
  bg = "white",
  res = 100)
image(log(Matrix.376), col=rev(COL2(diverging = "RdBu", n=200)), axes = FALSE)
dev.off()



