wd <- getwd()
homepath <- str_split_i(wd, "Normative_model", 1)
functionFolder.SCDev <- paste0(homepath, "/SC_development/Rcode_SCdevelopment/gamfunction")
source(paste0(functionFolder.SCDev, "/SCrankcorr.R"))
source(paste0(functionFolder.SCDev, "/gammsmooth.R"))
source(paste0(functionFolder.SCDev, '/gamderivatives.R'))
source(paste0(functionFolder.SCDev, "/plotdata_generate.R"))
compute_alignment <- function(dataname0, groupvar, TDlevel=1, ADHDlevel=0, ds.resolution){
  
  gam.data <- get(dataname0)
  
  # separate data by diagnosis
  sample_data_TD <<- gam.data[gam.data[[groupvar]]==TDlevel, ]
  sample_data_ADHD <<- gam.data[gam.data[[groupvar]]==ADHDlevel, ]
  
  # 1st. fitted developmental models
  #####################
  covariates<-"sex+meanFD"
  smooth_var<-"age"
  ## TD
  dataname<-"sample_data_TD"
  gamresult.TD <- list()
  gammod.TD <- list()
  for (i in 1:element_num){
    SClabel<-grep("SC.", names(sample_data_TD), value=T)[i]
    region<-SClabel
    gamresult<-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=FALSE, bootstrap=F)
    gamresult<-as.data.frame(gamresult)
    gamresult.TD[[i]] <- gamresult
    
    gammod <-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=T, bootstrap=F)
    gammod.TD[[i]] <- gammod
  }
  gamresult.TD.df <- do.call(rbind, lapply(gamresult.TD, function(x) data.frame(x)))
  gamresult.TD.df[,c(2:17)]<-lapply(gamresult.TD.df[,c(2:17)], as.numeric)
  
  ## ADHD
  dataname<-"sample_data_ADHD"
  gamresult.ADHD <- list()
  gammod.ADHD <- list()
  for (i in 1:element_num){
    SClabel<-grep("SC.", names(sample_data_ADHD), value=T)[i]
    region<-SClabel
    gamresult<-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=FALSE, bootstrap=F)
    gamresult<-as.data.frame(gamresult)
    gamresult.ADHD[[i]] <- gamresult
    
    gammod <-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=T, bootstrap=F)
    gammod.ADHD[[i]] <- gammod
  }
  gamresult.ADHD.df <- do.call(rbind, lapply(gamresult.ADHD, function(x) data.frame(x)))
  gamresult.ADHD.df[,c(2:17)]<-lapply(gamresult.ADHD.df[,c(2:17)], as.numeric)
  
  # alignment
  gamresult.TD.df <- gamresult.TD.df %>% mutate(
    meanderv2_2 = case_when(
      meanderv2 > mean(meanderv2)+3*sd(meanderv2) ~ NA,
      meanderv2 < mean(meanderv2)-3*sd(meanderv2) ~ NA,
      .default = meanderv2
    ))
  SCrank.TD <- SCrankcorr(gamresult.TD.df, "meanderv2_2", ds.resolution)
  SCrank.TD$diagnoselabel <- TDlevel
  
  gamresult.ADHD.df <- gamresult.ADHD.df %>% mutate(
    meanderv2_2 = case_when(
      meanderv2 > mean(meanderv2)+3*sd(meanderv2) ~ NA,
      meanderv2 < mean(meanderv2)-3*sd(meanderv2) ~ NA,
      .default = meanderv2
    ))
  SCrank.ADHD <- SCrankcorr(gamresult.ADHD.df, "meanderv2_2", ds.resolution)
  SCrank.ADHD$diagnoselabel <- ADHDlevel
  SCrank.df <- rbind(SCrank.TD, SCrank.ADHD)
  
  SCrank.TD.plotdf <- SCrankcorr(gamresult.TD.df, "meanderv2", ds.resolution, dsdata=T)
  SCrank.TD.plotdf$diagnoselabel <- TDlevel
  SCrank.ADHD.plotdf <- SCrankcorr(gamresult.ADHD.df, "meanderv2", ds.resolution, dsdata=T)
  SCrank.ADHD.plotdf$diagnoselabel <- ADHDlevel
  SCrank.plotdf <- rbind(SCrank.TD.plotdf, SCrank.ADHD.plotdf)
  
  # Other derivatives
  gamresults <- data.frame(parcel= rep(gamresult.TD.df$parcel, 2), increase.onset=c(gamresult.TD.df$increase.onset, 
                                                                                    gamresult.ADHD.df$increase.onset),
                           increase.offset=c(gamresult.TD.df$increase.offset, gamresult.ADHD.df$increase.offset),
                           diagnosis_label=c(rep(TDlevel, 78), rep(ADHDlevel, 78)))
  
  #####################
  
  # 2nd. compute 1st derivative and the alignment development trajectories.
  ## define SCrank
  SCrank <- matrix(NA, ds.resolution, ds.resolution)
  indexup.ds <- upper.tri(SCrank)
  indexsave.ds <- !indexup.ds
  for (x in 1:ds.resolution){
    for (y in 1:ds.resolution){
      SCrank[x,y] = x^2 + y^2
    }
  }
  SCrank <- rank(SCrank[indexsave.ds], ties.method = "average")
  
  ## compute derivatives
  derivative.sum.TD <- list()
  derivative.sum.ADHD <- list()
  draws<-1
  increments<-1000
  for (x in 1:element_num){
    # TD
    modobj<-gammod.TD[[x]]
    derivdata.TD <- gam.derivatives(modobj, "age", smoothvector= gam.data$age, draws, increments, return_posterior_derivatives = FALSE)
    derivdata.TD$label_ID<-gamresult.TD.df$parcel[x]
    plotdata.TD <- plotdata_generate(modobj, NA, "age", gam.data$age)
    firstfit <- plotdata.TD$fit[1]
    derivdata.TD$derivative.diw <- derivdata.TD$derivative / firstfit
    derivative.sum.TD[[x]] <- derivdata.TD
    
    # ADHD
    modobj<-gammod.ADHD[[x]]
    derivdata.ADHD <- gam.derivatives(modobj, "age", smoothvector= gam.data$age, draws, increments, return_posterior_derivatives = FALSE)
    derivdata.ADHD$label_ID<-gamresult.ADHD.df$parcel[x]
    plotdata.ADHD <- plotdata_generate(modobj, NA, "age", gam.data$age)
    firstfit <- plotdata.ADHD$fit[1]
    derivdata.ADHD$derivative.diw <- derivdata.ADHD$derivative / firstfit
    derivative.sum.ADHD[[x]] <- derivdata.ADHD
  }
  
  derivative.sum.TD.df <- do.call(rbind, lapply(derivative.sum.TD, function(x){t(x$derivative.diw)}))
  derivative.sum.ADHD.df <- do.call(rbind, lapply(derivative.sum.ADHD, function(x){t(x$derivative.diw)}))
  
  cor_coefs.TD <- apply(derivative.sum.TD.df, 2, function(col) cor(col, SCrank, method = "spearman"))
  cor_coefs.ADHD <- apply(derivative.sum.ADHD.df, 2, function(col) cor(col, SCrank, method = "spearman"))
  
  alignment.dev <- data.frame(age=derivdata.TD$age, cor_coefs.TD=as.numeric(t(cor_coefs.TD)), cor_coefs.ADHD=as.numeric(t(cor_coefs.ADHD)))
  
  resultlist <- list(SecDerivSCrank.df=SCrank.df, SecDerivSCrank.plotdf=SCrank.plotdf,gamresults=gamresults, alignment.dev=alignment.dev, 
                     derivative.sum.TD.df=derivative.sum.TD.df, derivative.sum.ADHD.df=derivative.sum.ADHD.df, SCrank=SCrank)
  
  return(resultlist)
}
