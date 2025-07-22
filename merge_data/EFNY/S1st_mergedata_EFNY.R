## This script is to generate a dataframe, in which each column is the strength for an edge.
## For schaefer 400 atlas, 70786 edges left after deleting edges connecting to limbic regions.
## The ID of edges which should be deleted according to consistency threshold will be saved out. 
library(R.matlab)
library(mgcv)
library(ggplot2)
library(tidyverse)
library(parallel)
library(reshape)
library(corrplot)
rm(list = ls())
wdpath <- getwd()
if (str_detect(wdpath, "cuizaixu_lab")){
  SC_path <-'/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/SCmat_Du15_all'
  demopath<-'/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/EFNY/demography'
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/EFNY/interfile'
  FigureFolder <- '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/EFNY/figures'
}else{
  # in PC
  demopath <- 'D:/xuxiaoyu/project_files/BP_youth_data/demography'
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/Normative_model/interfileFolder_EFNY'
  FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/Normative_model/Figures_EFNY'
  functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/Normative_model/functions'
}

# demo
demo_EFNY <- read.csv(paste0(demopath, '/basic_demo_merge_screen.csv')) # 152 subjects with complete dMRI & normal anat
log_file <- paste0(interfileFolder, "/missing_SC_data_log.txt")
write("", log_file)

# Du15 rank
Du15_rank <- read.csv('/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/atlas/Du15/Du15_SArank.csv')
Du15_rank <- Du15_rank %>%
  filter(NetworkLabel != "UNCERTAIN" & NetworkLabel != "NONE")
Du15_rank <- Du15_rank %>%
  mutate(rank = rank(MedianSA))
Du15_rank <- Du15_rank %>% 
  arrange(rank)


#################################################
#### import SC data
#### 15 regions, 16*15/2=120SCs
#################################################
colname <- character(length = 120)
for (i in 1:120){
  colname[i] <- paste0('SC.', as.character(i))
}
SCdata.sum<- data.frame(t(rep(0,120)))
names(SCdata.sum)<-colname
SCdata.sum$subID <- "NULL"


for (i in 1:nrow(demo_EFNY)){
  subID <- demo_EFNY$subID[i]
  SCname <- paste0(subID, '_sift_invnodevol_count.csv')
  if (file.exists(paste0(SC_path, '/', SCname))){
    SCmat <- read_csv(paste0(SC_path, '/', SCname), col_names = FALSE)
    SCmat <- as.matrix(SCmat)
    SCmat <- SCmat[-16, -16]
    
    SCmat <- SCmat[Du15_rank$X, Du15_rank$X]
    indexup <- upper.tri(SCmat)
    indexsave <- !indexup ###keep lower triangle and diagonal
    SCdat <- as.data.frame(c(SCmat[indexsave]))
    SCdat <- as.data.frame(t(SCdat), row.names = NULL)
    names(SCdat) <- colname
    row.names(SCdat) <- NULL
    SCdat$subID[1] <- subID
    SCdata.sum<-rbind(SCdata.sum, SCdat)
  }else{
    msg <- paste0('SC data for subject ', subID, ' does not exist.\n')
    cat(msg, file = log_file, append = TRUE)
  }
  
}

SCdata.sum<-SCdata.sum[-1,]
saveRDS(SCdata.sum, paste0(interfileFolder, '/EFNY.SCdata.sum.msmtcsd.delLM.rds'))
