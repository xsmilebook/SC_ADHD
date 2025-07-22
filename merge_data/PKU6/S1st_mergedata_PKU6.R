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
wdpath <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD"
if (str_detect(wdpath, "cuizaixu_lab")){
  SC_path <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/PKU6/SCmat_Du15_all'
  demopath<-'/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/PKU6/demography'
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/PKU6/interfile_PKU6'
  FigureFolder <- '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/PKU6/Figures_PKU6'
}else{
  # in PC
  SC_path <-'D:/code/R_projects/SC_ADHD/datasets/PKU6/SCmat_Du15_all'
  demopath <- 'D:/code/R_projects/SC_ADHD/datasets/PKU6/demography'
  interfileFolder <- 'D:/code/R_projects/SC_ADHD/datasets/PKU6/interfile_PKU6'
  FigureFolder<-'D:/code/R_projects/SC_ADHD/datasets/PKU6/Figures_PKU6'
  # functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/Normative_model/functions'
}

# demo
demo_PKU6 <- read.csv(paste0(demopath, '/basic_demo_merge.csv')) # 152 subjects with complete dMRI & normal anat


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


for (i in 1:nrow(demo_PKU6)){
  subID <- demo_PKU6$subID[i]
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
  }
}

SCdata.sum<-SCdata.sum[-1,]
saveRDS(SCdata.sum, paste0(interfileFolder, '/SCdata.sum.msmtcsd.delLM.rds'))
