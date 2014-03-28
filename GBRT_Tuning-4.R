##  UAH MS OR Capstone research
##  Kaggle AMS Solar Enger Prediction Competition
##
##  Jonathan Davis 2014-03-17
##
##  Script 4: Model tuning experiment



##load packages
library(ncdf4)
library(foreach)
library(reshape2)
library(doParallel)
library(leaps)
library(bestglm)
library(scales)
library(gbm)
library(ggplot2)
library(rgl)
library(caret)



#load workspace image
load("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")
setwd("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy")

setwd("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/")
load("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")


## get list of kriged data filenames
krigTrainfiles<-list.files("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTraining")
krigTrainfiles<-list.files("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTraining")



cl<-makeCluster(4) 
registerDoParallel(cl)

Mesonet<-foreach(i = 1:15, .combine = rbind, .packages = c("foreach", "ncdf4")) %dopar% {
  ##open connection to ncdf file
  GEFSfile<-nc_open(paste("Derived Data/KrigedTraining/", krigTrainfiles[i], sep = ""))
  GEFSvar<-ncvar_get(GEFSfile)
  Month<-GEFSfile$dim$Time$vals[1:5113]
  VarName<-gsub(".*Krig_(.*).nc*", "\\1", GEFSfile[[1]])
  nc_close(GEFSfile)
  
  vars<-foreach(l = 1:98, .combine = rbind) %do% {

    a1<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvar[l, 1, ], Var = paste0(VarName, 1))
    a2<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvar[l, 2, ], Var = paste0(VarName, 2))
    a3<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvar[l, 3, ], Var = paste0(VarName, 3))
    a4<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvar[l, 4, ], Var = paste0(VarName, 4))
    a5<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvar[l, 5, ], Var = paste0(VarName, 5))
#     
#     total<-data.frame(Time = 1:5113, Month = Month, Station = l, value = apply(data.frame(a1$value, a2$value, a3$value, a4$value, a5$value), 1, sum), Var = paste0(VarName, "TSum"))
#     Range<-data.frame(Time = 1:5113, Month = Month, Station = l, value = apply(data.frame(a1$value, a2$value, a3$value, a4$value, a5$value), 1, max) - apply(data.frame(a1$value, a2$value, a3$value, a4$value, a5$value), 1, min), Var = paste0(VarName, "TRange"))
#     
    
    rbind(a1, a2, a3, a4, a5)
  }
  vars
  
}

stopCluster(cl)

##correct for negative values that can result from the Kriging estimates in some cases for accumulated precip
Mesonet$value[Mesonet$value < 0]<-0

MesonetC<-dcast(Mesonet, Station + Time + Month ~ Var )
rm(Mesonet)
gc()
MesonetC$Month<-factor(floor(((MesonetC$Month)%%8766)/730.5) +1)


## REplace station id's with Station names
MesonetC$Station<-stationinfo$stid[MesonetC$Station]

## melt the training data to get Y values mergable with Mesonet
train$Date<-1:5113

trainMelt<-melt(id = "Date", measured = names(train)[2:99], data = train)
MesonetC<-merge(MesonetC, trainMelt, by.x = c("Time", "Station"), by.y = c("Date", "variable"))
rm(trainMelt)
gc()

##  create formula object
#fullform<-as.formula(paste(paste("value ~ (", paste0(names(Mesonet)[6:75], collapse = " + "), ")^2 +", sep = ""), paste0("I(", names(Mesonet)[6:75], "^2)", collapse = " + "), sep = ""))
form<-as.formula(paste("value ~ ", paste0(sample(names(MesonetC)[5:(dim(MesonetC)[2] -1)][regexpr(".*[1-5]$", names(MesonetC)[5:(dim(MesonetC)[2] -1)])>0], length(names(MesonetC)[5:(dim(MesonetC)[2] -1)][regexpr(".*[1-5]$", names(MesonetC)[5:(dim(MesonetC)[2] -1)])>0])), collapse = " + ")))
formwAgg<-as.formula(paste("value ~ ", paste0(sample(names(MesonetC)[5:(dim(MesonetC)[2] -1)], length(names(MesonetC)[5:(dim(MesonetC)[2] -1)])), collapse = " + ")))

rm(Mesonet)
gc()
 
# #residual full model
# fullform<-as.formula(paste(paste("Res ~ (", paste0(names(Mesonet)[6:75], collapse = " + "), ")^2 +", sep = ""), paste0("I(", names(Mesonet)[6:75], "^2)", collapse = " + "), sep = ""))
# additive<-as.formula(paste("Res ~ ", paste0(names(Mesonet)[6:75], collapse = " + ")))


##  get list of station names
stationnames<-stationinfo$stid
S<-sample(1:98, 1)
i = S

stationnamesR<-stationinfo$stid[S]

set.seed(13747)
##  store a seed for each ncdf4 file
seeds<-floor(2^10*runif(98))

##  set script seed with a purely random number regurgitated straight from my brain

idepth<-c(5,10,15)

learn<-c(.05, .03, .01)
ntree<-c(300, 800, 2000)

cvfolds<-c(2,3,6)

gbm.Models<-foreach(cv = 1:3, .combine = rbind) %:% foreach(lr = 1:3, .combine = rbind) %:% 
  foreach(d = 1:3, .combine = rbind) %:% foreach(i = 1:length(S), .combine = rbind, .packages = c("gbm")) %do% {

  print(paste(cv, lr, d))  
  set.seed(seeds[i + 10])

  ##  extract dats for station i 
  dat1<-MesonetC[MesonetC[,"Station"] == stationnames[i],]
  dat1<-dat1[order(dat1[,"Time"]), ]
  dat<-dat1[1:3800,]
  dat1<-dat1[3801:5113,]
  dat$r<-runif(dim(dat)[1])
  dat<-dat[order(dat[,"r"]),]

   
    ## create boosted regression trees
  rt<-system.time({
  GBR1<-gbm(form, data = dat, bag.fraction = .5, distribution = "laplace", train.fraction = .67, cv.folds = cvfolds[cv],interaction.depth =idepth[d], shrinkage = learn[lr], n.trees = ntree[lr],  n.cores = 3)
  best.iter1<-gbm.perf(GBR1, method = "cv", plot.it = T)
  predgbm1<-predict.gbm(GBR1, newdata = dat1, n.trees = best.iter1)
  })
  rm(GBR1)
  gc()
  
  rt<-as.numeric(rt[3])
  
  GBR2<-gbm(form, data = dat, bag.fraction = .5, distribution = 'laplace', train.fraction = .67, cv.folds = cvfolds[cv],interaction.depth = idepth[d], shrinkage = learn[lr], n.trees = ntree[lr],  n.cores = 3)
  best.iter2<-gbm.perf(GBR2, method = "cv", plot.it = F)
  predgbm2<-predict.gbm(GBR2, newdata = dat1, n.trees = best.iter2)
  rm(GBR2)
  gc()
  
  GBR3<-gbm(form, data = dat, bag.fraction = .5, distribution = 'laplace', train.fraction = .67, cv.folds = cvfolds[cv],interaction.depth = idepth[d], shrinkage = learn[lr], n.trees = ntree[lr],  n.cores = 3)
  best.iter3<-gbm.perf(GBR3, method = "cv", plot.it = F)
  predgbm3<-predict.gbm(GBR3, newdata = dat1, n.trees = best.iter3)
  rm(GBR3)
  gc()
  
  GBR4<-gbm(form, data = dat, bag.fraction = .5, distribution = 'laplace', train.fraction = .67, cv.folds = cvfolds[cv],interaction.depth = idepth[d], shrinkage = learn[lr], n.trees = ntree[lr],  n.cores = 3)
  best.iter4<-gbm.perf(GBR4, method = "cv", plot.it = F)
  predgbm4<-predict.gbm(GBR4, newdata = dat1, n.trees = best.iter4)
  rm(GBR4)
  gc()
  
  GBR5<-gbm(form, data = dat, bag.fraction = .5, distribution = 'laplace', train.fraction = .67, cv.folds = cvfolds[cv],interaction.depth = idepth[d], shrinkage = learn[lr], n.trees = ntree[lr],  n.cores = 3)
  best.iter5<-gbm.perf(GBR5, method = "cv", plot.it = F)
  predgbm5<-predict.gbm(GBR5, newdata = dat1, n.trees = best.iter5)
  rm(GBR5)
  gc()
  
  
  GBR6<-gbm(form, data = dat, bag.fraction = .5, distribution = 'laplace', train.fraction = .67, cv.folds = cvfolds[cv],interaction.depth = idepth[d], shrinkage = learn[lr], n.trees = ntree[lr],  n.cores = 3)
  best.iter6<-gbm.perf(GBR6, method = "cv", plot.it = F)
  predgbm6<-predict.gbm(GBR6, newdata = dat1, n.trees = best.iter6)
  rm(GBR6)
  gc()
  
  
  GBR7<-gbm(form, data = dat, bag.fraction = .5, distribution = 'laplace', train.fraction = .67, cv.folds = cvfolds[cv],interaction.depth = idepth[d], shrinkage = learn[lr], n.trees = ntree[lr],  n.cores = 3)
  best.iter7<-gbm.perf(GBR7, method = "cv", plot.it = F)
  predgbm7<-predict.gbm(GBR7, newdata = dat1, n.trees = best.iter6)
  rm(GBR7)
  gc()
  
  GBR8<-gbm(form, data = dat, bag.fraction = .5, distribution = 'laplace', train.fraction = .67, cv.folds = cvfolds[cv],interaction.depth = idepth[d], shrinkage = learn[lr], n.trees = ntree[lr],  n.cores = 3)
  best.iter8<-gbm.perf(GBR8, method = "cv", plot.it = T)
  predgbm8<-predict.gbm(GBR8, newdata = dat1, n.trees = best.iter6)
  rm(GBR8)
  gc()
  
  
  
  
    pred<-predgbm1
    pred2<-apply(data.frame(predgbm1, predgbm2), 1, mean)
    pred3<-apply(data.frame(predgbm1, predgbm2, predgbm3), 1, mean)  
    pred4<-apply(data.frame(predgbm1, predgbm2, predgbm3, predgbm4), 1, mean)  
    pred5<-apply(data.frame(predgbm1, predgbm2, predgbm3, predgbm4, predgbm5), 1, mean)  
    pred6<-apply(data.frame(predgbm1, predgbm2, predgbm3, predgbm4, predgbm5, predgbm6), 1, mean)  
    pred7<-apply(data.frame(predgbm1, predgbm2, predgbm3, predgbm4, predgbm5, predgbm6, predgbm7), 1, mean)  
    pred8<-apply(data.frame(predgbm1, predgbm2, predgbm3, predgbm4, predgbm5, predgbm6, predgbm7, predgbm8), 1, mean)  
  
    GBRMAE<-mean(abs(pred - dat1$value))
    GBRMAE2<-mean(abs(pred2 - dat1$value))
    GBRMAE3<-mean(abs(pred3 - dat1$value))
    GBRMAE4<-mean(abs(pred4 - dat1$value))
    GBRMAE5<-mean(abs(pred5 - dat1$value))
    GBRMAE6<-mean(abs(pred6 - dat1$value))
    GBRMAE7<-mean(abs(pred7 - dat1$value))
    GBRMAE8<-mean(abs(pred8 - dat1$value))
  
  
    a1<-data.frame(Station = stationnames[i], Loss = GBRMAE, Aggregated = 1, CVFolds = cvfolds[cv], InteractionDepth = idepth[d], shrinkage = learn[lr], runtime = rt)                                                                          
    a2<-data.frame(Station = stationnames[i], Loss = GBRMAE2, Aggregated = 2, CVFolds = cvfolds[cv], InteractionDepth = idepth[d], shrinkage = learn[lr], runtime = 2*rt)                                                                          
    a3<-data.frame(Station = stationnames[i], Loss = GBRMAE3, Aggregated = 3, CVFolds = cvfolds[cv], InteractionDepth = idepth[d], shrinkage = learn[lr], runtime = 3*rt)                                                                          
    a4<-data.frame(Station = stationnames[i], Loss = GBRMAE4, Aggregated = 4, CVFolds = cvfolds[cv], InteractionDepth = idepth[d], shrinkage = learn[lr], runtime = 4*rt)                                                                          
    a5<-data.frame(Station = stationnames[i], Loss = GBRMAE5, Aggregated = 5, CVFolds = cvfolds[cv], InteractionDepth = idepth[d], shrinkage = learn[lr], runtime = 5*rt)                                                                          
    a6<-data.frame(Station = stationnames[i], Loss = GBRMAE6, Aggregated = 6, CVFolds = cvfolds[cv], InteractionDepth = idepth[d], shrinkage = learn[lr], runtime = 6*rt)                                                                          
    a7<-data.frame(Station = stationnames[i], Loss = GBRMAE7, Aggregated = 7, CVFolds = cvfolds[cv], InteractionDepth = idepth[d], shrinkage = learn[lr], runtime = 7*rt)                                                                          
    a8<-data.frame(Station = stationnames[i], Loss = GBRMAE8, Aggregated = 8, CVFolds = cvfolds[cv], InteractionDepth = idepth[d], shrinkage = learn[lr], runtime = 8*rt)                                                                          

    res<-rbind(a1, a2, a3, a4, a5, a6, a7, a8)
    res                       
 


} ##end i 

save.image("./Derived Data/Figures/GBMTuning3.RData")

