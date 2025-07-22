## This script is to generate a dataframe, in which each column is the strength for an edge in large-scale network.
## For schaefer 400 --> Yeo network atlas, elementnum mean edges in large-scale network.
library(R.matlab)
library(ggplot2)
library(tidyverse)
library(parallel)
library(openxlsx)
library(rjson)
library(reshape)
library(psych)
rm(list = ls())
homepath <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD"
# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
elementnum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2

SC_path_EFNY <-'/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/defaultatlas'
SC_path_PKU6 <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/PKU6/SCmat'
SC_path_CCNP <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/SC'
Volume_path_EFNY <-'/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/schaefer400_nodevolume'
Volume_path_PKU6 <-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/PKU6/schaefer400_nodevolume'
Volume_path_CCNP <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/CCNP/processed/schaefer400_nodevolume'

demopath <- paste0(homepath, '/datasets/demography')
interfileFolder <- paste0(homepath, '/datasets/interfileFolder')
functionFolder <- paste0(homepath, "/code/functions")
resultFolder <- paste0(homepath, "/datasets/results")
FigureFolder <- paste0(homepath, "/Figures/Yeo17")

Behavior <- read.csv(paste0(demopath, '/THREE_SITES_basic_demo.csv'))
# Behavior <- Behavior %>% distinct(ID, .keep_all = TRUE)

#### import schaefer400 index
schaefer400_index_SA<-read.csv(paste0(homepath, '/datasets/atlas/schaefer400/schaefer400_index_SA.csv'))
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

# filter index of P75th and P25th of CV.
deleteindex75 <- readRDS(paste0(interfileFolder, '/CV75_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))
deleteindex25 <- readRDS(paste0(interfileFolder, '/CV25_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))

# assign each region to Yeo.resolution*Yeo.resolution network.
schaefer400_index.Yeo7 <- schaefer400_index[order(schaefer400_index$index_7network_LRmixed),]
schaefer400_index.Yeo17 <- schaefer400_index[order(schaefer400_index$index_17network_LRmixed),]
schaefer400_index.Yeo7 <- schaefer400_index.Yeo7 %>% mutate(Yeo.resolutionnode = recode_factor(network_label,
                                                                                               "Vis" = 1,
                                                                                               "SomMot" = 2,
                                                                                               "DorsAttn" = 3,
                                                                                               "SalVentAttn" = 4,
                                                                                               "Cont" = 5,
                                                                                               "Default" = 6,
                                                                                               "Limbic" = 0))

print(summary(schaefer400_index.Yeo7$Yeo.resolutionnode))
schaefer400_index.Yeo17 <- schaefer400_index.Yeo17 %>% 
  mutate(Yeo.resolutionnode = recode_factor(network_label_17network,
                                            "VisCent" = 3,
                                            "VisPeri" = 1,
                                            "SomMotA" =2,
                                            "SomMotB" = 4,
                                            "DorsAttnA" =5,
                                            "DorsAttnB" = 6,
                                            "SalVentAttnA" =9,
                                            "SalVentAttnB" =12,
                                            "ContA" =11,
                                            "ContB" =14,
                                            "ContC" =7,
                                            "DefaultA" = 13,
                                            "DefaultB" = 15,
                                            "DefaultC" = 8,
                                            "TempPar" = 10))

print(summary(schaefer400_index.Yeo17$Yeo.resolutionnode))

if (Yeoresolution == 7){
  Yeo.resolutionnode <- schaefer400_index.Yeo7$Yeo.resolutionnode
  Yeo.resolutionnode <- factor(Yeo.resolutionnode, levels = c(1:Yeoresolution.delLM))
}else if (Yeoresolution == 17){
  Yeo.resolutionnode <- schaefer400_index.Yeo17$Yeo.resolutionnode
  Yeo.resolutionnode <- factor(Yeo.resolutionnode, levels = c(1:Yeoresolution.delLM))
}else{
  print("Invalid Yeoresolution!")
}

# SC 376*376 --> Yeoresolution.delLM*Yeoresolution.delLM
matrixYeo.resolution <- matrix(NA, Yeoresolution.delLM, Yeoresolution.delLM)
matrixYeo.resolution[lower.tri(matrixYeo.resolution, diag = T)] <- c(1:elementnum)
matrixYeo.resolution[upper.tri(matrixYeo.resolution)] <- t(matrixYeo.resolution)[upper.tri(matrixYeo.resolution)]
matrix_SCYeo.resolution <- matrix(NA, 376, 376)
for (x in 1:Yeoresolution.delLM){
  for (y in 1:Yeoresolution.delLM){
    xindex <- which(Yeo.resolutionnode==x)
    yindex <- which(Yeo.resolutionnode==y)
    matrix_SCYeo.resolution[xindex, yindex] <- matrixYeo.resolution[x,y]
  }
}
# an index telling how 376*376 map to 12*12
Yeo.resolution.index <- matrix_SCYeo.resolution[lower.tri(matrix_SCYeo.resolution, diag = T)]
#################################################

#### import SC data
#### Yeoresolution.delLM regions, (Yeoresolution.delLM+1)*Yeoresolution.delLM/2=elementnum SCs
#### extract a dataframe containing elementnum columns, each represents an edge.
#################################################
colname <- character(length = elementnum)
for (i in 1:elementnum){
  colname[i] <- paste0('SC.', as.character(i))
}

if (str_detect(homepath, "cuizaixu_lab")){
  SCdata.sum <- mclapply(1:nrow(Behavior), function(i){
    result <- tryCatch({
      scanID <- Behavior$scanID[i]
      site <- Behavior$site[i]
      if (site == "EFNY"){
        SCname <- paste0(scanID, '_dir-PA_space-T1w_desc-preproc_msmtconnectome.mat')
        SC_file_path <- paste0(SC_path_EFNY, '/', SCname)
        volumefile <- paste0(Volume_path_EFNY, '/', scanID, '_Volume7.txt')
      }else if (site == "PKU6"){
        SCname <- paste0(scanID, '.mat')
        SC_file_path <- paste0(SC_path_PKU6, '/', SCname)
        volumefile <- paste0(Volume_path_PKU6, '/', scanID, '_Volume7.txt')
      }else{
        SCname <- paste0(scanID, '_space-T1w_desc-preproc_dhollanderconnectome.mat')
        SC_file_path <- file.path(SC_path_CCNP, SCname)
        volumefile <- paste0(Volume_path_CCNP, '/', scanID, '_Volume7.txt')
      }
      
      
      if (file.exists(SC_file_path)){
        SCmat <- readMat(SC_file_path)
        # load steamline counts matrix & fiber length matrix
        SCmat_raw <- SCmat$schaefer400.sift.radius2.count.connectivity[schaefer376_delLM, schaefer376_delLM]
        length_raw <- SCmat$schaefer400.radius2.meanlength.connectivity[schaefer376_delLM, schaefer376_delLM]
        if (Yeoresolution == 7){
          SCmat_raw <- SCmat_raw[orderYeo_7, orderYeo_7] # 376*376 nodes sorted by Yeo index
          length_raw <- length_raw[orderYeo_7, orderYeo_7] # 376*376 nodes sorted by Yeo index
        }else if (Yeoresolution == 17){
          SCmat_raw <- SCmat_raw[orderYeo_17, orderYeo_17]
          length_raw <- length_raw[orderYeo_17, orderYeo_17]
        }else{
          print("Invalid Yeoresolution!")
        }
        
        totallength_raw <- length_raw * SCmat_raw
        indexup <- upper.tri(SCmat_raw)
        indexsave <- !indexup
        SCmat_raw <- SCmat_raw[indexsave] # 1*70876 each element represents streamline counts
        SCmat_raw75 <- SCmat_raw25 <- SCmat_raw
        SCmat_raw75[deleteindex75]<-0 # remove top 1/4 inconsistent connetions
        SCmat_raw25[deleteindex25]<-0 # remove top 3/4 inconsistent connetions
        totallength_raw <- totallength_raw[indexsave]
        totallength_raw75 <- totallength_raw25 <- totallength_raw
        totallength_raw75[deleteindex75]<-0
        totallength_raw25[deleteindex25]<-0
        df <- data.frame(
          group = Yeo.resolution.index,
          value75 = SCmat_raw75,
          value25 = SCmat_raw25,
          length75 = totallength_raw75,
          length25 = totallength_raw25
        )
        # compute the sum of streamline counts / length for each fraction, in total of elementnum.
        result <- df %>%
          group_by(group) %>%
          summarise(sum_value75 = sum(value75), sum_value25 = sum(value25), sum_length75=sum(length75), 
                    sum_length25 = sum(length25))
        mean_length75 <- (result$sum_length75 / result$sum_value75)[1:elementnum]
        mean_length25 <- (result$sum_length25 / result$sum_value25)[1:elementnum]
        sumSC.raw75 <- result$sum_value75[1:elementnum]
        sumSC.raw25 <- result$sum_value25[1:elementnum]
        ## node volume
        # MARK 1: 
        if (any(is.na(sumSC.raw75))) {
          scanID <- Behavior$scanID[i]
          print(paste("!!! WARNING for scanID:", scanID, "-> sumSC.raw75 contains NA values!"))
        }
        if (file.exists(volumefile)){
          nodevolume <- read_table(volumefile, col_names=F, n_max = 400)
          if (nrow(nodevolume)==400){
            nodevolume <- as.numeric(nodevolume$X2[schaefer376_delLM]) # delete limbic regions

            if (Yeoresolution == 7){
              nodevolume <- nodevolume[orderYeo_7] # sorted by Yeo index
            }else if (Yeoresolution == 17){
              nodevolume <- nodevolume[orderYeo_17]
            }else{
              print("Invalid Yeoresolution!")
            }
            
            df2 <- data.frame(
              group = Yeo.resolutionnode,
              value = nodevolume
            )
            result2 <- df2 %>% arrange(group) %>%
              group_by(group) %>% 
              summarise(sum_value = sum(value))
            nodevolume_sum <- result2$sum_value[1:Yeoresolution.delLM] # sum of nodes' volume for each node fraction (Yeo.resolution).
            
            # MARK 3:
            if (any(is.na(nodevolume_sum)) || any(nodevolume_sum == 0)) {
              scanID <- Behavior$scanID[i]
              print(paste("!!! WARNING for scanID:", scanID, "-> nodevolume_sum contains NA or zero!"))
              print("nodevolume_sum values:")
              print(nodevolume_sum)
            }

            ### Yeo.resolution*Yeo.resolution
            volumeSC <- matrix(NA, Yeoresolution.delLM, Yeoresolution.delLM)
            for (x in 1:Yeoresolution.delLM){
              for (y in 1:Yeoresolution.delLM){
                volumeSC[x,y] <- (nodevolume_sum[x]+nodevolume_sum[y])/2
              }
            }
            volumeSC <- volumeSC[lower.tri(volumeSC, diag = T)] # the scale values of node volume for each edge.
            sumSC.invnode75 <- sumSC.raw75 / volumeSC
            sumSC.invnode25 <- sumSC.raw25 / volumeSC
          }else{
            sumSC.invnode75 <- rep(NA, elementnum)
            sumSC.invnode25 <- rep(NA, elementnum)
          }}else{
            sumSC.invnode75 <- rep(NA, elementnum)
            sumSC.invnode25 <- rep(NA, elementnum)
          }
        ###keep lower triangle and diagonal
        SCdat75 <- as.data.frame(sumSC.invnode75)
        SCdat75 <- as.data.frame(t(SCdat75), row.names = NULL)
        names(SCdat75) <- colname
        row.names(SCdat75) <- NULL
        SCdat75$scanID[1] <- scanID
        #SCdata.sum75<-rbind(SCdata.sum75, SCdat75)
        
        SCdat25 <- as.data.frame(sumSC.invnode25)
        SCdat25 <- as.data.frame(t(SCdat25), row.names = NULL)
        names(SCdat25) <- colname
        row.names(SCdat25) <- NULL
        SCdat25$scanID[1] <- scanID
        #SCdata.sum25<-rbind(SCdata.sum25, SCdat25)
        
        # not inverse node volume
        SCdat75_noInvNode <- as.data.frame(sumSC.raw75)
        SCdat75_noInvNode <- as.data.frame(t(SCdat75_noInvNode), row.names = NULL)
        names(SCdat75_noInvNode) <- colname
        row.names(SCdat75_noInvNode) <- NULL
        SCdat75_noInvNode$scanID[1] <- scanID
        
        SCdat25_noInvNode <- as.data.frame(sumSC.raw25)
        SCdat25_noInvNode <- as.data.frame(t(SCdat25_noInvNode), row.names = NULL)
        names(SCdat25_noInvNode) <- colname
        row.names(SCdat25_noInvNode) <- NULL
        SCdat25_noInvNode$scanID[1] <- scanID
        # MARK 4:
        if (all(is.na(sumSC.invnode75))) {
          scanID <- Behavior$scanID[i]
          print(paste(">>> FINAL CHECK for scanID:", scanID, "-> All sumSC.invnode75 are NA <<<"))
        }
      }else{
        print(paste("File not found for scanID:", scanID, "at path:", SC_file_path))
      }
      print(i)
      return(data.frame(SCdat75, SCdat25, SCdat75_noInvNode, SCdat25_noInvNode))
    }, error = function(e){
      scanID <- Behavior$scanID[i] # 尝试获取ID
      print(paste("### ERROR processing scanID:", scanID, "at iteration", i, "###"))
      print(paste("Error message:", e$message))
      return(NULL) # 返回 NULL，表示这个被试处理失败
    })
    return(result)
  },  mc.cores = 50)
  print("The parallel Mat reading is over!")
  saveRDS(SCdata.sum, paste0(interfileFolder, "/SCdata.sum.list_Yeo", Yeoresolution, ".rds"))
  
  SCdata.sum75 <- do.call(rbind, lapply(SCdata.sum, function(x) as.data.frame(x[,c(1:(elementnum+1))])))
  SCdata.sum25 <- do.call(rbind, lapply(SCdata.sum, function(x) as.data.frame(x[,c((elementnum+2):(elementnum*2+2))])))
  names(SCdata.sum25) <- c(colname, "scanID")
  
  SCdata.sum75_noInvNode <- do.call(rbind, lapply(SCdata.sum, function(x) as.data.frame(x[,c((elementnum*2+3):(elementnum*3+3))])))
  SCdata.sum25_noInvNode <- do.call(rbind, lapply(SCdata.sum, function(x) as.data.frame(x[,c((elementnum*3+4):(elementnum*4+4))])))
  names(SCdata.sum25_noInvNode) <- c(colname, "scanID")
  names(SCdata.sum75_noInvNode) <- c(colname, "scanID")
  
  SCdata.sum75.merge <- merge(SCdata.sum75, Behavior, by="scanID")
  SCdata.sum25.merge <- merge(SCdata.sum25, Behavior, by="scanID")
  SCdata.sum75_noInvNode.merge <- merge(SCdata.sum75_noInvNode, Behavior, by="scanID")
  SCdata.sum25_noInvNode.merge <- merge(SCdata.sum25_noInvNode, Behavior, by="scanID")
  
  # convert variables' types
  SCdata.sum75.merge$scanID <- as.factor(SCdata.sum75.merge$scanID) ; SCdata.sum75.merge$site <- as.factor(SCdata.sum75.merge$site)
  SCdata.sum25.merge$scanID <- as.factor(SCdata.sum25.merge$scanID) ; SCdata.sum25.merge$site <- as.factor(SCdata.sum25.merge$site)
  SCdata.sum75_noInvNode.merge$scanID <- as.factor(SCdata.sum75_noInvNode.merge$scanID) ; SCdata.sum75_noInvNode.merge$site <- as.factor(SCdata.sum75_noInvNode.merge$site)
  SCdata.sum25_noInvNode.merge$scanID <- as.factor(SCdata.sum25_noInvNode.merge$scanID) ; SCdata.sum25_noInvNode.merge$site <- as.factor(SCdata.sum25_noInvNode.merge$site)
  
  saveRDS(SCdata.sum75.merge, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
  saveRDS(SCdata.sum25.merge, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV25_sumSCinvnode.sum.msmtcsd.merge.rds'))
  saveRDS(SCdata.sum75_noInvNode.merge, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSC.sum.msmtcsd.merge.rds'))
  saveRDS(SCdata.sum25_noInvNode.merge, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV25_sumSC.sum.msmtcsd.merge.rds'))
}else{
  SCdata.sum75.merge <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge.rds'))
  SCdata.sum25.merge <- readRDS(paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV25_sumSCinvnode.sum.msmtcsd.merge.rds'))
}

# plot
Matrix.tmp <- matrix(NA, nrow = Yeoresolution.delLM, ncol=Yeoresolution.delLM)
linerange_frame<-data.frame(x=c(0.5,Yeoresolution.delLM+0.5), ymin =rep(-Yeoresolution.delLM-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -Yeoresolution.delLM-0.5), xmin=rep(0.5, times=2), xmax=rep(Yeoresolution.delLM+0.5, times=2))
SCdata.sum75.merge <- SCdata.sum75.merge[!is.na(SCdata.sum75.merge$SC.19),]
SCdata.sum25.merge <- SCdata.sum25.merge[!is.na(SCdata.sum25.merge$SC.19),]
# 75
age.time <- c(9, 11, 13)
meanSCdata.sepage75 <- list()
for (i in 1:3){
  age.tmp <- age.time[i]
  index.tmp <- which(SCdata.sum75.merge$age>=age.tmp & SCdata.sum75.merge$age<(age.tmp+1))
  SCdata.tmp <- SCdata.sum75.merge[index.tmp, ]
  meanSCdata.tmp <- colMeans(SCdata.tmp[,c(2:(elementnum+1))])
  meanSCdata.sepage75[[i]] <- meanSCdata.tmp}

SCmin <- min(c(meanSCdata.sepage75[[1]], meanSCdata.sepage75[[2]], meanSCdata.sepage75[[3]]))
SCmax <- max(c(meanSCdata.sepage75[[1]], meanSCdata.sepage75[[2]], meanSCdata.sepage75[[3]]))
for (i in 1:3){
  age.tmp <- age.time[i]
  index.tmp <- which(SCdata.sum75.merge$age>=age.tmp & SCdata.sum75.merge$age<(age.tmp+1))
  SCdata.tmp <- SCdata.sum75.merge[index.tmp, ]
  meanSCdata.tmp <- colMeans(SCdata.tmp[,c(2:(elementnum+1))])
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- meanSCdata.tmp
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, Yeoresolution.delLM)
  rownames(Matrix.tmp) <-seq(1, Yeoresolution.delLM)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, Yeoresolution.delLM)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)

  Fig<-ggplot(data =matrixtmp.df.melt)+
    geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
    scale_fill_distiller(type="seq", palette = "RdBu",limit=c(SCmin-0.1, SCmax+0.1), na.value = "grey")+
    scale_color_distiller(type="seq", palette = "RdBu",limit=c(SCmin-0.1, SCmax+0.1),na.value = "grey")+
    geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
    geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
    geom_segment(aes(x = 0.5 , y = -0.5 , xend = Yeoresolution.delLM+0.5 ,yend = -Yeoresolution.delLM-0.5), color="black", linewidth=0.5)+
    ggtitle(label = paste("CV75, age", age.tmp, "~", (age.tmp+1)))+labs(x=NULL, y=NULL)+
    scale_y_continuous(breaks=NULL, labels = NULL)+
    scale_x_continuous(breaks=NULL, labels = NULL)+
    theme(axis.line = element_blank(),
          #axis.ticks=element_line(linewidth = 0),
          axis.text.x=element_text(size=Yeoresolution.delLM, angle=45, hjust=1),
          axis.text.y=element_text(size=Yeoresolution.delLM, angle=315, hjust=1,vjust=1),
          axis.title =element_text(size=18),
          plot.title = element_text(size=18, hjust = 0.5),
          legend.title=element_text(size=18),
          legend.text=element_text(size=18),
          panel.background=element_rect(fill=NA),
          panel.grid.major=element_line(linewidth = 0),
          panel.grid.minor=element_line(linewidth = 1))
  Fig
  filename<-paste0(FigureFolder,"/CV75/Matrix_Yeo",Yeoresolution, "_sumSCinvnode_Age8_22/Age_", age.tmp, "_ds.resolutionnet_delLM_CV75.tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm", create.dir = T)

}

# 25
age.time <- c(9, 11, 13)
meanSCdata.sepage25 <- list()
for (i in 1:3){
  age.tmp <- age.time[i]
  index.tmp <- which(SCdata.sum25.merge$age>=age.tmp & SCdata.sum25.merge$age<(age.tmp+1))
  SCdata.tmp <- SCdata.sum25.merge[index.tmp, ]
  meanSCdata.tmp <- colMeans(SCdata.tmp[,c(2:(elementnum+1))])
  meanSCdata.sepage25[[i]] <- meanSCdata.tmp}

SCmin <- min(c(meanSCdata.sepage25[[1]], meanSCdata.sepage25[[2]], meanSCdata.sepage25[[3]]))
SCmax <- max(c(meanSCdata.sepage25[[1]], meanSCdata.sepage25[[2]], meanSCdata.sepage25[[3]]))

for (i in 1:3){
  age.tmp <- age.time[i]
  index.tmp <- which(SCdata.sum25.merge$age>=age.tmp & SCdata.sum25.merge$age<(age.tmp+1))
  SCdata.tmp <- SCdata.sum25.merge[index.tmp, ]
  meanSCdata.tmp <- colMeans(SCdata.tmp[,c(2:(elementnum+1))])
  meanSCdata.sepage25[[i]] <- meanSCdata.tmp
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- meanSCdata.tmp
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, Yeoresolution.delLM)
  rownames(Matrix.tmp) <-seq(1, Yeoresolution.delLM)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, Yeoresolution.delLM)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)

  Fig<-ggplot(data =matrixtmp.df.melt)+
    geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
    scale_fill_distiller(type="seq", palette = "RdBu",limit=c(SCmin-0.1, SCmax+0.1), na.value = "grey")+
    scale_color_distiller(type="seq", palette = "RdBu",limit=c(SCmin-0.1, SCmax+0.1),na.value = "grey")+
    geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
    geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
    geom_segment(aes(x = 0.5 , y = -0.5 , xend = Yeoresolution.delLM+0.5 ,yend = -Yeoresolution.delLM-0.5), color="black", linewidth=0.5)+
    ggtitle(label = paste("CV25, age", age.tmp, "~", (age.tmp+1)))+labs(x=NULL, y=NULL)+
    scale_y_continuous(breaks=NULL, labels = NULL)+
    scale_x_continuous(breaks=NULL, labels = NULL)+
    theme(axis.line = element_blank(),
          #axis.ticks=element_line(linewidth = 0),
          axis.text.x=element_text(size=Yeoresolution.delLM, angle=45, hjust=1),
          axis.text.y=element_text(size=Yeoresolution.delLM, angle=315, hjust=1,vjust=1),
          axis.title =element_text(size=18),
          plot.title = element_text(size=18, hjust = 0.5),
          legend.title=element_text(size=18),
          legend.text=element_text(size=18),
          panel.background=element_rect(fill=NA),
          panel.grid.major=element_line(linewidth = 0),
          panel.grid.minor=element_line(linewidth = 1))
  Fig
  filename<-paste0(FigureFolder, "/CV25/Matrix_Yeo",Yeoresolution, "_sumSCinvnode_Age8_22/Age_", age.tmp, "_",Yeoresolution.delLM, "net_delLM_CV25.tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm", create.dir = T)

}



