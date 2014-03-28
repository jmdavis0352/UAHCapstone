##  UAH MS OR Capstone research
##  Kaggle AMS Solar Enger Prediction Competition
##
##  Jonathan Davis 2014-03-20
##
##  Script 5: Calculate preditcions over the test period



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



#load workspace image
load("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")
load("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")
load("c:/users/jmdavis_2/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")

setwd("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/")
setwd("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/")
setwd("c:/users/jmdavis_2/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/")

## get list of kriged data filenames
krigTrainfiles<-list.files("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTrainingSE")
krigTrainfiles<-list.files("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTrainingSE")

krigTestfiles<-list.files("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTest5SE")
krigTrainfiles<-list.files("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTrainingSE")

# ## get list of kriged data filenames
# krigTrainfiles<-list.files("c:/users/jmdavis_2/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTraining")
# krigTestfiles<-list.files("c:/users/jmdavis_2/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTest5")



############  Compile the training data ##################################

cl<-makeCluster(4) 
registerDoParallel(cl)

Mesonet<-foreach(i = 1:15, .combine = rbind, .packages = c("foreach", "ncdf4")) %dopar% {
  ##open connection to ncdf file
  GEFSfile<-nc_open(paste("Derived Data/KrigedTrainingSE/", krigTrainfiles[i], sep = ""),)
  GEFSvar<-ncvar_get(GEFSfile,GEFSfile$var[[1]]$name)
  GEFSvarSE<-ncvar_get(GEFSfile,GEFSfile$var[[2]]$name)
  Month<-GEFSfile$dim$Time$vals[1:5113]
  VarName<-gsub(".*Krig_(.*).nc*", "\\1", GEFSfile[[1]])
  VarNameSE<-paste0(gsub(".*Krig_(.*).nc*", "\\1", GEFSfile[[1]]), "SE")
  nc_close(GEFSfile)
  
  vars<-foreach(l = 1:98, .combine = rbind) %do% {

    a1<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvar[l, 1, ], Var = paste0(VarName, 1))
    a2<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvar[l, 2, ], Var = paste0(VarName, 2))
    a3<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvar[l, 3, ], Var = paste0(VarName, 3))
    a4<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvar[l, 4, ], Var = paste0(VarName, 4))
    a5<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvar[l, 5, ], Var = paste0(VarName, 5))
#

    rbind(a1, a2, a3, a4, a5 )
  }


#   varsSE<-foreach(l = 1:98, .combine = rbind) %do% {
#   
#   a1<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvarSE[l, 1, ], Var = paste0(VarName, 1, "SE"))
#   a2<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvarSE[l, 2, ], Var = paste0(VarName, 2, "SE"))
#   a3<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvarSE[l, 3, ], Var = paste0(VarName, 3, "SE"))
#   a4<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvarSE[l, 4, ], Var = paste0(VarName, 4, "SE"))
#   a5<-data.frame(Time = 1:5113, Month = Month, Station =l,  value = GEFSvarSE[l, 5, ], Var = paste0(VarName, 5, "SE"))
#   #
#   
#   rbind(a3, a4, a5 )
#   }
# 
# rbind(vars, varsSE)


vars

  
}

stopCluster(cl)

##correct for negative values that can result from the Kriging estimates in some cases for accumulated precip
Mesonet$value[Mesonet$value < 0]<-0

MesonetC<-dcast(Mesonet, Station + Time + Month ~ Var )
rm(Mesonet)
gc()
MesonetC$Month<-factor(floor(((MesonetC$Month)%%8766)/730.5) +1)


## Add station info
MesonetC$Lat<-stationinfo$nlat[MesonetC$Station]
MesonetC$Lon<-stationinfo$elon[MesonetC$Station]
MesonetC$Elev<-stationinfo$elev[MesonetC$Station]
MesonetC$Station<-stationinfo$stid[MesonetC$Station]


## melt the training data to get Y values mergable with Mesonet
train$Date<-1:5113

trainMelt<-melt(id = "Date", measured = names(train)[2:99], data = train)
MesonetC<-merge(MesonetC, trainMelt, by.x = c("Time", "Station"), by.y = c("Date", "variable"))
rm(trainMelt)
gc()

#############  Compile the test data ####################################

cl<-makeCluster(4) 
registerDoParallel(cl)

MesonetTest<-foreach(i = 1:15, .combine = rbind, .packages = c("foreach", "ncdf4")) %dopar% {
  ##open connection to ncdf file
  GEFSfile<-nc_open(paste("Derived Data/KrigedTest5SE/", krigTestfiles[i], sep = ""))
  GEFSvar<-ncvar_get(GEFSfile, GEFSfile$var[[1]]$name)
  GEFSvarSE<-ncvar_get(GEFSfile,GEFSfile$var[[2]]$name)
  Month<-GEFSfile$dim$Time$vals[1:1796]
  VarName<-gsub(".*Krig_(.*).nc*", "\\1", GEFSfile[[1]])
  nc_close(GEFSfile)
  
  vars<-foreach(l = 1:98, .combine = rbind) %do% {
    
    a1<-data.frame(Time = 1:1796, Month = Month, Station =l,  value = GEFSvar[l, 1, ], Var = paste0(VarName, 1))
    a2<-data.frame(Time = 1:1796, Month = Month, Station =l,  value = GEFSvar[l, 2, ], Var = paste0(VarName, 2))
    a3<-data.frame(Time = 1:1796, Month = Month, Station =l,  value = GEFSvar[l, 3, ], Var = paste0(VarName, 3))
    a4<-data.frame(Time = 1:1796, Month = Month, Station =l,  value = GEFSvar[l, 4, ], Var = paste0(VarName, 4))
    a5<-data.frame(Time = 1:1796, Month = Month, Station =l,  value = GEFSvar[l, 5, ], Var = paste0(VarName, 5))
    
     
    rbind(a1, a2, a3, a4, a5)
  }
  
#   varsSE<-foreach(l = 1:98, .combine = rbind) %do% {
#     
#     a1<-data.frame(Time = 1:1796, Month = Month, Station =l,  value = GEFSvarSE[l, 1, ], Var = paste0(VarName, 1, "SE"))
#     a2<-data.frame(Time = 1:1796, Month = Month, Station =l,  value = GEFSvarSE[l, 2, ], Var = paste0(VarName, 2, "SE"))
#     a3<-data.frame(Time = 1:1796, Month = Month, Station =l,  value = GEFSvarSE[l, 3, ], Var = paste0(VarName, 3, "SE"))
#     a4<-data.frame(Time = 1:1796, Month = Month, Station =l,  value = GEFSvarSE[l, 4, ], Var = paste0(VarName, 4, "SE"))
#     a5<-data.frame(Time = 1:1796, Month = Month, Station =l,  value = GEFSvarSE[l, 5, ], Var = paste0(VarName, 5, "SE"))
#     #
#     
#     rbind(a3, a4, a5 )
#   }
  
#  rbind(vars, varsSE)
  
vars

}

stopCluster(cl)

##correct for negative values that can result from the Kriging estimates in some cases for accumulated precip
MesonetTest$value[MesonetTest$value<0]<-0

MesonetCT<-dcast(MesonetTest, Station + Time + Month ~ Var )
rm(MesonetTest)
gc()
MesonetCT$Month<-factor(floor(((MesonetCT$Month)%%8766)/730.5) +1)


## REplace station id's with Station names
MesonetCT$Lat<-stationinfo$nlat[MesonetCT$Station]
MesonetCT$Lon<-stationinfo$elon[MesonetCT$Station]
MesonetCT$Elev<-stationinfo$elev[MesonetCT$Station]
MesonetCT$Station<-stationinfo$stid[MesonetCT$Station]




##  create formula object
#fullform<-as.formula(paste(paste("value ~ (", paste0(names(Mesonet)[6:75], collapse = " + "), ")^2 +", sep = ""), paste0("I(", names(Mesonet)[6:75], "^2)", collapse = " + "), sep = ""))
form<-as.formula(paste("value ~ ", paste0(sample(names(MesonetC)[3:(dim(MesonetC)[2] -1)][regexpr(".*SE$", names(MesonetC)[3:(dim(MesonetC)[2] -1)])<0], length(names(MesonetC)[3:(dim(MesonetC)[2] -1)][regexpr(".*SE$", names(MesonetC)[3:(dim(MesonetC)[2] -1)])<0])), collapse = " + ")))
formSE<-as.formula(paste("value ~ ", paste0(sample(names(MesonetC)[3:(dim(MesonetC)[2] -1)], length(names(MesonetC)[3:(dim(MesonetC)[2] -1)])), collapse = " + ")))

MesonetC$r<-runif(dim(MesonetC)[1])
MesonetC<-MesonetC[order(MesonetC[,"r"]),]
# 
# #residual full model
# fullform<-as.formula(paste(paste("Res ~ (", paste0(names(Mesonet)[6:75], collapse = " + "), ")^2 +", sep = ""), paste0("I(", names(Mesonet)[6:75], "^2)", collapse = " + "), sep = ""))
# additive<-as.formula(paste("Res ~ ", paste0(names(Mesonet)[6:75], collapse = " + ")))


##  get list of station names
stationnames<-stationinfo$stid


set.seed(13747)
##  store a seed for each ncdf4 file
seeds<-floor(2^10*runif(98))
##  read the dates from the submission file so i don't have to copy and paste them over and over
submissionDate<-read.csv("Source Data/submissionDate.csv", header = TRUE)

predlist<-list(98)

allModel<-gbm(form, data = MesonetC, bag.fraction = .5, distribution = "laplace", train.fraction = .67, cv.folds = 3,interaction.depth = 15, shrinkage = .1, n.trees = 900,  n.cores = 3)
optTrees<-gbm.perf(allModel, method = "cv", plot.it = T)
allModel$cv.error[optTrees]

allModel1<-gbm(form, data = MesonetC, bag.fraction = .5, distribution = "laplace", train.fraction = .67, cv.folds = 3,interaction.depth = 15, shrinkage = .1, n.trees = 900,  n.cores = 3)
optTrees1<-gbm.perf(allModel1, method = "cv", plot.it = T)
allModel1$cv.error[optTrees1]

save.image("TwoModels600.Rdata")

allModel2<-gbm(form, data = MesonetC, bag.fraction = .5, distribution = "laplace", train.fraction = .67, cv.folds = 3,interaction.depth = 15, shrinkage = .1, n.trees = 900,  n.cores = 3)
optTrees2<-gbm.perf(allModel2, method = "cv", plot.it = T)
allModel2$cv.error[optTrees2]

allModel3<-gbm(form, data = MesonetC, bag.fraction = .5, distribution = "laplace", train.fraction = .67, cv.folds = 3,interaction.depth = 15, shrinkage = .1, n.trees = 900,  n.cores = 3)
optTrees3<-gbm.perf(allModel3, method = "cv", plot.it = T)
allModel3$cv.error[optTrees3]


gbm.Models<-foreach(i = 1:98, .combine = rbind, .packages = c("gbm"), .errorhandling = 'stop') %do% {

  print(paste(i))  


test<-MesonetCT[MesonetCT[,"Station"] == stationnames[i],]
  test<-test[order(test[,"Time"]), ]
#    
  pred1<-predict.gbm(allModel,newdata = test, n.trees = optTrees)
  pred2<-predict.gbm(allModel1,newdata = test, n.trees = optTrees1)
  pred3<-predict.gbm(allModel2,newdata = test, n.trees = optTrees2)


     pred<-apply(data.frame(pred1, pred2, pred3), 1, mean)
#     fitted<-apply(data.frame(fitted1, fitted2, fitted3, fitted4, fitted5, fitted6, fitted7), 1, mean)
#    
#     a<-list(submission = data.frame(Station = stationnames[i], Date = 1:1796, value = pred, stringsAsFactors = FALSE),
#                               fitted = data.frame(Station = stationnames[i], Time = 1:5113, value = fitted)  )
       
    submission <- data.frame(Station = stationnames[i], Date = 1:1796, value = pred, stringsAsFactors = FALSE)
    submission




} ##end i 



## cast VPredict to produce a submission file formatted matrix
vPredictCast<-dcast(gbm.Models, Date ~ Station, value.var = 'value')
vPredictCast$Date<-submissionDate$Date


write.csv(vPredictCast, "./Derived Data/Submission.csv", row.names = FALSE)

