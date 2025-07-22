library(R.matlab)
library(psych)
source("D:/xuxiaoyu/GeneralRfunctions/rotate_parcellation/perm.sphere.p.R")
source("D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction/permspheregam.R")

SCrankcorr <- function(gamresult, computevar, ds.resolution, perm.id.full, dsdata=FALSE,
                       setgam=FALSE, set_fx=FALSE){

  #### connectional rank
  ########################
  #connectional rank
  Matrix.ds<-matrix(NA, nrow=ds.resolution, ncol=ds.resolution)
  indexup.ds <- upper.tri(Matrix.ds)
  indexsave.ds <- !indexup.ds
  Matrix.ds.index<-Matrix.ds
  SClength.ds=ds.resolution*(ds.resolution+1)/2 # number of elements in the lower triangle of mat
  Matrix.ds.index[indexsave.ds]<-c(1:SClength.ds) # initial order of index in lower triangle
  Matrix.ds.SCrank<-Matrix.ds
  for (x in 1:ds.resolution){
    for (y in 1:ds.resolution){
      Matrix.ds.SCrank[x,y]<-x^2+y^2
    }
  }
  Matrix.ds.SCrank[indexup.ds]<-NA
  Matrix.ds.SCrank[indexsave.ds]<-rank(Matrix.ds.SCrank[indexsave.ds], ties.method = "average")
  gamresult.ds<-data.frame(SCrank=Matrix.ds.SCrank[indexsave.ds], computevar=NA)
  gamresult.ds$computevar <-gamresult[, computevar]
  
  correstimate<-corr.test(gamresult.ds, method="spearman")$r[2,1]
  p.spearman<-corr.test(gamresult.ds, method="spearman")$p[2,1]
  # compute null rho using perm ID from left and right hemisphere separately.
  perm.id.ds.L <- perm.id.full[1:ds.resolution, ]
  perm.id.ds.R <- perm.id.full[(ds.resolution+1):(ds.resolution*2), ]-ds.resolution
  perm.id.ds <- cbind(perm.id.ds.L, perm.id.ds.R)
  # get perm ID for lower triangle of matrix
  SC.perm.id.ds<-matrix(data=NA, nrow = SClength.ds, ncol=20000)
  SC.perm.id.sep <- list()
  for (i in 1:20000){
    ordertmp<-perm.id.ds[,i]
    tmp<-Matrix.ds.index[ordertmp, ordertmp] # permutate nodes and index of lower triangle
    for (x in 1:ds.resolution){
      for (y in 1:ds.resolution){
        if (is.na(tmp[x,y])){
          tmp[x,y]<-tmp[y,x]}
      }}
    SC.perm.id.sep[[i]] <- tmp[indexsave.ds]
  }# return permutated index of lower triangle
  
  for (j in 1:20000){
    SC.perm.id.ds[,j]<-SC.perm.id.sep[[j]] # cbind perm ID for lower triangle of matrix
  }
  ## if there is NA in gamresult.ds$computevar, then remove the corresponding row in perm ID.
  if (length(which(is.na(gamresult.ds$computevar)))>0){
    SC.perm.id.ds <- SC.perm.id.ds[-which(is.na(gamresult.ds$computevar)),]
    SC.perm.id.ds.test<-lapply(c(1:20000), function(x) rank(SC.perm.id.ds[,x], ties.method="first"))
    for (i in 1:20000){
      SC.perm.id.ds[,i]<-SC.perm.id.ds.test[[i]]
    }
    gamresult.ds.noNA <- gamresult.ds[-which(is.na(gamresult.ds$computevar)), ]
  }else{gamresult.ds.noNA <- gamresult.ds}
  
  if (setgam==FALSE){
    pspin <- perm.sphere.p(gamresult.ds.noNA$computevar, gamresult.ds.noNA$SCrank, SC.perm.id.ds,corr.type='spearman')
  }else{
    pspin <- perm.sphere.p.gam(gamresult.ds.noNA$SCrank, gamresult.ds.noNA$computevar, SC.perm.id.ds, set_fx)
  }
  
  SCrankdf <- data.frame(ds.resolution=ds.resolution, Interest.var=computevar,
                         r.spearman=correstimate, p.spin=pspin, p.spearman=p.spearman)
  if (dsdata == TRUE){
    names(gamresult.ds) <- c("SCrank", computevar)
    return(gamresult.ds)
  }else{
    return(SCrankdf)
  }
  
}



