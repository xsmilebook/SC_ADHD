library(mgcv)
perm.sphere.p.gam = function(x,y,perm.id, setfx=FALSE) {
  
  nroi = dim(perm.id)[1]  # number of regions
  nperm = dim(perm.id)[2] # number of permutations
  
  gammodel.emp <- gam(y ~ s(x, k=3, fx=setfx), method="REML")
  gamresult.emp <- summary(gammodel.emp)
  Rsq.emp = gamresult.emp[["r.sq"]]  # empirical correlation
  
  # permutation of measures
  x.perm = y.perm = array(NA,dim=c(nroi,nperm))
  for (r in 1:nperm) {
    for (i in 1:nroi) {
      x.perm[i,r] = x[perm.id[i,r]]
      y.perm[i,r] = y[perm.id[i,r]]
    }
  }
  
  # correlation to unpermuted measures
  Rsq.null.xy = Rsq.null.yx = vector(length=nperm)
  for (r in 1:nperm) {
    x.tmp <- x.perm[,r]
    gammodel.x <- gam(y ~ s(x.tmp, k=3, fx=setfx), method="REML")
    gamresult.x <- summary(gammodel.x)
    Rsq.null.xy[r] <- gamresult.x[["r.sq"]]
    
    y.tmp <- y.perm[,r]
    gammodel.y <- gam(y.tmp ~ s(x, k=3, fx=setfx), method="REML")
    gamresult.y <- summary(gammodel.y)
    Rsq.null.yx[r] <- gamresult.y[["r.sq"]]
  }
  
  # p-value definition depends on the sign of the empirical correlation
  if (Rsq.emp>0) {
    p.perm.xy = sum(Rsq.null.xy>Rsq.emp)/nperm
    p.perm.yx = sum(Rsq.null.yx>Rsq.emp)/nperm
  }
  
  # return average p-value
  return((p.perm.xy+p.perm.yx)/2)
  
}