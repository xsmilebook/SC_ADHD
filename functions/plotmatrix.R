# plot RdBu ds.resolution*ds.resolution matrix plot
library(reshape)
library(ggplot2)

plotmatrix <- function(dataname, variable, ds.resolution, Pvar=NA, NAcol="white", lmthr=NA, 
                       axeslabels=NULL, axeslabelsGap=T, linerange_frame=NA, PaletteSet=NA, Pvar.noFDR=NA){
  data.tmp <- get(dataname)
  Matrix.tmp <- matrix(NA, ds.resolution, ds.resolution)
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- data.tmp[[variable]]
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, ds.resolution)
  rownames(Matrix.tmp) <-seq(1, ds.resolution)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, ds.resolution)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  
  if (is.na(Pvar)==0){
    if (sum(data.tmp[[Pvar]]<0.05, na.rm=T) > 0){
      Matrix.tmp.sig <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
      Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- (data.tmp[[Pvar]]<0.05)
      Matrix.tmp.sig[upper.tri(Matrix.tmp.sig)] <- t(Matrix.tmp.sig)[upper.tri(Matrix.tmp.sig)]
      colnames(Matrix.tmp.sig) <-seq(1, ds.resolution)
      rownames(Matrix.tmp.sig) <-seq(1, ds.resolution)
      matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
      matrixtmp.df.sig$nodeid <- seq(1, ds.resolution)
      matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig,id.vars=c("nodeid"))
      matrixtmp.df.sig.melt$variable<-as.numeric(matrixtmp.df.sig.melt$variable)
      matrixtmp.df.sig.melt$nodeid<-0-matrixtmp.df.sig.melt$nodeid
      matrixtmp.df.sig.melt$value<-as.numeric(matrixtmp.df.sig.melt$value)
      matrixtmp.df.sig.melt <- matrixtmp.df.sig.melt[-which(matrixtmp.df.sig.melt$value==0),]
      textcolor <- "black"
    }else{
      Matrix.tmp.sig <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
      matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
      matrixtmp.df.sig$nodeid <- seq(1, ds.resolution)
      matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig,id.vars=c("nodeid"))
      matrixtmp.df.sig.melt$variable<-as.numeric(matrixtmp.df.sig.melt$variable)
      matrixtmp.df.sig.melt$nodeid<-0-matrixtmp.df.sig.melt$nodeid
      matrixtmp.df.sig.melt$value<-as.numeric(matrixtmp.df.sig.melt$value)
      textcolor <- "transparent"
    }
  }else{
    Matrix.tmp.sig <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
    matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
    matrixtmp.df.sig$nodeid <- seq(1, ds.resolution)
    matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig,id.vars=c("nodeid"))
    matrixtmp.df.sig.melt$variable<-as.numeric(matrixtmp.df.sig.melt$variable)
    matrixtmp.df.sig.melt$nodeid<-0-matrixtmp.df.sig.melt$nodeid
    matrixtmp.df.sig.melt$value<-as.numeric(matrixtmp.df.sig.melt$value)
    textcolor <- "transparent"
  }
  
  if (is.na(Pvar.noFDR)==0){
    # add significance not survive FDR correction.
    Matrix.tmp.sig_noFDR <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
    Matrix.tmp.sig_noFDR[lower.tri(Matrix.tmp.sig_noFDR, diag = T)] <- (data.tmp[[Pvar.noFDR]] <0.05 & data.tmp[[Pvar]] > 0.05)
    Matrix.tmp.sig_noFDR[upper.tri(Matrix.tmp.sig_noFDR)] <- t(Matrix.tmp.sig_noFDR)[upper.tri(Matrix.tmp.sig_noFDR)]
    colnames(Matrix.tmp.sig_noFDR) <-seq(1, ds.resolution)
    rownames(Matrix.tmp.sig_noFDR) <-seq(1, ds.resolution)
    matrixtmp.df.sig_noFDR <- as.data.frame(Matrix.tmp.sig_noFDR)
    matrixtmp.df.sig_noFDR$nodeid <- seq(1, ds.resolution)
    matrixtmp.df.sig_noFDR.melt <- melt(matrixtmp.df.sig_noFDR,id.vars=c("nodeid"))
    matrixtmp.df.sig_noFDR.melt$variable<-as.numeric(matrixtmp.df.sig_noFDR.melt$variable)
    matrixtmp.df.sig_noFDR.melt$nodeid<-0-matrixtmp.df.sig_noFDR.melt$nodeid
    matrixtmp.df.sig_noFDR.melt$value<-as.numeric(matrixtmp.df.sig_noFDR.melt$value)
    matrixtmp.df.sig_noFDR.melt <- matrixtmp.df.sig_noFDR.melt[-which(matrixtmp.df.sig_noFDR.melt$value==0),]
  }else{
    Matrix.tmp.sig_noFDR <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
    matrixtmp.df.sig_noFDR <- as.data.frame(Matrix.tmp.sig_noFDR)
    matrixtmp.df.sig_noFDR$nodeid <- seq(1, ds.resolution)
    matrixtmp.df.sig_noFDR.melt <- melt(matrixtmp.df.sig_noFDR,id.vars=c("nodeid"))
    matrixtmp.df.sig_noFDR.melt$variable<-as.numeric(matrixtmp.df.sig_noFDR.melt$variable)
    matrixtmp.df.sig_noFDR.melt$nodeid<-0-matrixtmp.df.sig_noFDR.melt$nodeid
    matrixtmp.df.sig_noFDR.melt$value<-as.numeric(matrixtmp.df.sig_noFDR.melt$value)
  }
  
  
  plottitle <- paste(dataname, variable)
  if (is.na(lmthr)){
    lmthr <- max(abs(data.tmp[[variable]]))
  }
  
  if (length(linerange_frame)==1){
    linerange_frame<-data.frame(x=c(0.5,ds.resolution+0.5), ymin =rep(-ds.resolution-0.5, times=2), ymax =rep(-0.5, times=2),
                                y=c(-0.5, -ds.resolution-0.5), xmin=rep(0.5, times=2), xmax=rep(ds.resolution+0.5, times=2))
  }
  
    
  if (is.null(axeslabels)){
    axesbreaks <- NULL
    axesbreaksy <- NULL
  }else if (axeslabelsGap==F){
    axesbreaks <- c(1:ds.resolution)
    axesbreaksy <- -1*c(1:ds.resolution)
  }else{
    axesbreaks <- seq(from=2, to=ds.resolution, by=2)
    if (ds.resolution %% 2 == 0){axisyend = (ds.resolution-1)}else{axisyend = ds.resolution}
    axesbreaksy <- seq(from=-1, to=-axisyend, by=-2)
  }
  
  if (length(PaletteSet) == 1){
    PaletteSet <- list(Name="RdBu", drirection=-1, lmmin = -lmthr, lmmax = lmthr)
  }
  
  Fig <- ggplot(data=matrixtmp.df.melt)+
    geom_tile(aes(x=variable, y=nodeid, fill=value, color=value))+
    scale_fill_distiller(type="seq", palette = PaletteSet$Name, limit=c(PaletteSet$lmmin, PaletteSet$lmmax),direction=PaletteSet$drirection, na.value = NAcol)+
    scale_color_distiller(type="seq", palette = PaletteSet$Name,limit=c(PaletteSet$lmmin, PaletteSet$lmmax),direction=PaletteSet$drirection, na.value = NAcol)+
    geom_text(data =matrixtmp.df.sig_noFDR.melt, aes(x=variable, y=nodeid, label = "+"), color=textcolor, vjust = 0.5, hjust = 0.5, size=8)+
    geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), color=textcolor, vjust = 0.8, hjust = 0.5, size=10)+
    geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
    geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
    geom_segment(aes(x = 0.5 , y = -0.5 , xend = ds.resolution+0.5 ,yend = -ds.resolution-0.5), color="black", linewidth=0.5)+
    ggtitle(label = plottitle)+labs(x=NULL, y=NULL)+
    scale_y_continuous(breaks=axesbreaksy, labels = axeslabels[abs(axesbreaksy)])+
    scale_x_continuous(breaks=axesbreaks, labels = axeslabels[abs(axesbreaks)], position = "bottom")+
    theme(axis.line = element_blank(), 
          #axis.ticks=element_line(linewidth = 0),
          axis.text.x=element_text(size=27, angle=90, hjust=1, vjust=0.5, color="black"), 
          axis.text.y=element_text(size=27, angle=0, hjust=1,vjust=0.5, color="black"),
          axis.title =element_text(size=18),
          plot.title = element_text(size=18, hjust = 0.5),
          legend.title=element_text(size=18),
          legend.text=element_text(size=18), aspect.ratio = 1,
          panel.background=element_rect(fill=NA),
          panel.grid.major=element_line(linewidth = 0), 
          panel.grid.minor=element_line(linewidth = 1))
  
  return(Fig)
}
