library(R.matlab)
library(psych)
source("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/perm.sphere.p.R")
source("/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction/permspheregam.R")

SCrankcorr <- function(gamresult, matsize, computevar, ds.resolution, perm.id.full, rhonull=FALSE,
                       set_fx=FALSE){
  Matrix<-matrix(NA, nrow=matsize, ncol=matsize)
  indexup<-upper.tri(Matrix)
  indexsave<-!indexup
  
  #### connectional rank
  ########################
  #connectional rank
  Matrix.ds<-matrix(NA, nrow=ds.resolution, ncol=ds.resolution)
  indexup.ds <- upper.tri(Matrix.ds)
  indexsave.ds <- !indexup.ds
  Matrix.ds.index<-Matrix.ds
  SClength.ds=ds.resolution*(ds.resolution+1)/2
  Matrix.ds.index[indexsave.ds]<-c(1:SClength.ds)
  Matrix.ds.SCrank<-Matrix.ds
  for (x in 1:ds.resolution){
    for (y in 1:ds.resolution){
      Matrix.ds.SCrank[x,y]<-(x+y)^2+(x-y)^2
    }
  }
  Matrix.ds.SCrank[indexup.ds]<-NA
  Matrix.ds.SCrank[indexsave.ds]<-rank(Matrix.ds.SCrank[indexsave.ds], ties.method = "average")
  gamresult.ds<-data.frame(SCrank=Matrix.ds.SCrank[indexsave.ds], computevar=NA)
  gamresult.ds$computevar <-gamresult[, computevar]
  
  correstimate<-corr.test(gamresult.ds, method="spearman")$r[2,1]
  
  ## spin test
  perm.id.ds <- matrix(NA, nrow=ds.resolution, ncol=10000)
  for (i in 1:ds.resolution){
    x0=floor((i-1)*matsize/ds.resolution)+1
    xt=floor((i)*matsize/ds.resolution)
    if (x0==xt){
      perm.id.ds[i, ]<- perm.id.full[x0, ]
    }else{
      perm.id.ds[i, ]<- colMeans(perm.id.full[x0:xt, ])
    }
  }
  
  perm.id.ds.tmp<-lapply(c(1:10000), function(x) rank(perm.id.ds[,x], ties.method="first"))
  for (i in 1:10000){
    perm.id.ds[,i]<-perm.id.ds.tmp[[i]]
  }
  SC.perm.id.ds<-matrix(data=NA, nrow = SClength.ds, ncol=10000)
  SC.perm.id.sep<-mclapply(1:10000, function(i){
    ordertmp<-perm.id.ds[,i]
    tmp<-Matrix.ds.index[ordertmp, ordertmp]
    for (x in 1:ds.resolution){
      for (y in 1:ds.resolution){
        if (is.na(tmp[x,y])){
          tmp[x,y]<-tmp[y,x]}
      }}
    return(tmp[indexsave.ds])}, mc.cores = 4)
  for (j in 1:10000){
    SC.perm.id.ds[,j]<-SC.perm.id.sep[[j]]
  }
  ## if NA
  if (length(which(is.na(gamresult.ds$computevar)))>0){
    SC.perm.id.ds <- SC.perm.id.ds[-which(is.na(gamresult.ds$computevar)),]
    SC.perm.id.ds.test<-lapply(c(1:10000), function(x) rank(SC.perm.id.ds[,x], ties.method="first"))
    for (i in 1:10000){
      SC.perm.id.ds[,i]<-SC.perm.id.ds.test[[i]]
    }
    gamresult.ds.noNA <- gamresult.ds[-which(is.na(gamresult.ds$computevar)), ]
  }else{gamresult.ds.noNA <- gamresult.ds}
  
  pspin <- perm.sphere.p(gamresult.ds.noNA$computevar, gamresult.ds.noNA$SCrank, SC.perm.id.ds,corr.type='spearman')
  
  rhonullvalue <- perm.sphere.p(gamresult.ds.noNA$computevar, gamresult.ds.noNA$SCrank, SC.perm.id.ds,corr.type='spearman', rhonull=TRUE)
  SCrankdf <- data.frame(ds.resolution=ds.resolution, matsize=matsize, Interest.var=computevar,
                         r.spearman=correstimate, p.spin=pspin)
  
  if (rhonull == TRUE){
    return(rhonullvalue)
  }else{
    return(SCrankdf)
  }
  
}



