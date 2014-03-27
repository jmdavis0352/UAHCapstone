##  UAH MS OR Capstone research
##  Kaggle AMS Solar Enger Prediction Competition
##
##  Jonathan Davis 2014-02-25
##
##  Script 1: reads csv files and extracts some basic information for use in subsequent scripts. 

##
##  the ncdf4 library was not available through CRAN on windows when i did this anlaysis, so you have to download the zip file from the authors website:
##  http://cirrus.ucsd.edu/~pierce/ncdf/
##  there are detailed instructions for installation. On linux you'll need to have Netcdf libraries already installed before you try to install ncdf4 or you'll get a compile error. This is pretty 
##  easy. Go to you favorite package manager and sert netcdf and several obvious choices will come up. I just kept clicking them until it worked. 


##  load libraries
library(ncdf4)

##set working directory to an appropriate location on your machine
setwd("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/")

##  read csv files adjust file paths accoring to where you store the source data
stationinfo<-read.csv("Source Data/station_info.csv", header = TRUE, stringsAsFactors = FALSE)
train<-read.csv("Source Data/train_dependent.csv", header = TRUE, stringsAsFactors = FALSE)


##  get list of file names from directory where source data is stored 
trainfiles<-list.files("Source Data/train_independent/")
testfiles<-list.files("Source Data/test_independent/")

##  open one file to extract list of lat long values
GEFSfile<-nc_open(paste("Source Data/train_independent/", trainfiles[1], sep = ""))
Lats<-as.numeric(GEFSfile$dim$lat$vals)
Longs<-as.numeric(GEFSfile$dim$lon$vals)
nc_close(GEFSfile)
rm(GEFSfile)

##  write out workspace image
save.image("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")

