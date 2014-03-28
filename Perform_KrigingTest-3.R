##  UAH MS OR Capstone research
##  Kaggle AMES Solar Enger Prediction Competition
##
##  Jonathan Davis 2014-02-25
##
##  Script 2: Perform Krigings on GEFs data to the Mesonet station locations. Write this data back into the NCDF files for 
##  efficient extraction during model development scripts 

##  load libraries
library(lhs)
library(ncdf4)
library(geoR)
library(scales)
library(IDPmisc)
library(rgl)
library(foreach)
library(ggplot2)
library(geospt)
library(doParallel)
library(Rsolnp)

##  set working directory
setwd("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/")
setwd("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/")
setwd("c:/users/jmdavis_2/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/")


load("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")
load("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")
load("c:/users/jmdavis_2/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")





##  set script seed with a purely random number regurgitated straight from my brain
set.seed(13747)
##  store a seed for each ncdf4 file
seeds<-floor(2^12*runif(length(trainfiles)))
##  store locations of Mesonet stations for Kriging 
stationlocs<-stationinfo[ c("nlat", "elon")]

########################## Run parallel ########################## 
cl<-makeCluster(2) 
registerDoParallel(cl)

##  Begin loop over all 15 files 
Kriged<-foreach(f = 1:15, .combine = rbind, .packages = c("lhs", "geoR", "ncdf4", "foreach"),
                .errorhandling = 'stop' ) %dopar% {
  ##set random seed for this file
  set.seed(seeds[f])
  ##  Extract ncdf4 variable data
  GEFSfile<-nc_open(paste("Source Data/test_independent/", testfiles[f], sep = ""))
  GEFSvar<-ncvar_get(GEFSfile)
  VarName<-gsub(".*ent/(.*)_lat.*", "\\1", GEFSfile[[1]])
  ##  Create a data frame for storing the mean value 
  ncdf4data<-data.frame(lat = numeric(144), lon = numeric(144), Value = numeric(144))
  
  ##  Create vectors of indices over which to loop and extract data for kriging these should be the same
  ##  for each file, but just to be sure i'll reinitialize them on each one (only happens 15 times)
  Lats<-as.numeric(GEFSfile$dim$lat$vals)
  Longs<-as.numeric(GEFSfile$dim$lon$vals)
  Ensemble<-as.numeric(GEFSfile$dim$ens$vals)
  Hours<-as.numeric(GEFSfile$dim$fhour$vals)
  Times<-as.numeric(GEFSfile$dim$time$vals)
    
  ##  close connection to ncdf4 file
  nc_close(GEFSfile)
  
  ## create new ncdf4 file for storing interpolated GEFS Data
  dimStation<-ncdim_def(name ="Station", units ="Stations", vals = 1:98)
  dimHour<-ncdim_def(name ="Hour", units ="hours", vals = c(1:5))
  dimTime<-ncdim_def(name ="Time", units ="days", vals = Times)
  
  ##  create ncdf4 variable 
  krigGEFSVAR<-ncvar_def(name =VarName, units = "", dim = list(dimStation, dimHour,  dimTime), prec = 'double')
  krigGEFSVARSE<-ncvar_def(name =paste0(VarName, "SE"), units = "", dim = list(dimStation, dimHour,  dimTime), prec = 'double')
  
  ##  create file for storing this variable.
  #krigGEFfile<-nc_create(paste0("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTraining/train_Krig_",VarName, ".nc"), krigGEFSVAR)
  #krigGEFfile<-nc_create(paste0("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTest/test_Krig_",VarName, ".nc"), krigGEFSVAR)
  krigGEFfile<-nc_create(paste0("c:/users/jmdavis_2/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/KrigedTest5SE/test_Krig_",VarName, ".nc"), list(krigGEFSVAR,krigGEFSVARSE) )
  
  ##  loop over time of day 3-5, and times 1:5113 to create kriging predictions at the Mesonet locations length(Times)
  Results<-foreach(d= 1:5, .combine = rbind, .errorhandling = 'remove') %:% foreach(t = 1:length(Times), .combine = rbind, .verbose = FALSE, .errorhandling = 'remove') %do% {
    print(paste(f,d, t))
    ##R  andomly select points to remove from kriging fit near the center of the grid (area of interest)
    lhs<-optimumLHS(10, 2,eps = .05 ) ##generate design
    lhs[,1]<-ceiling(lhs[,1]/(1/10)) + 256 ##convert to discrete longitude values
    lhs[,2]<-ceiling(lhs[,2]/(1/4)) + 33 ##convert to discrete lattitude values
    lhs<-unique(data.frame(lon = lhs[,1], lat = lhs[,2], stringsAsFactors = FALSE)) ##unique because the rounding can cause duplicates
    lhs$Include<-"No"   
    
    ##  extract GEFS data
    for (lat in seq_along(Lats))  {
      for (long in seq_along(Longs)){
        val<-sum(GEFSvar[long, lat, d, , t])/11 ## get ensemble mean at long, lat, d, t
        ncdf4data[(lat-1)*length(Longs) + long, ] <-c(Lats[lat], Longs[long], val )
                
      }##end long
    }##end lat
    
    ## merge the ncdf4 data with lhs to extract the train and test data points
    krigingGrid<-merge(ncdf4data, lhs, by.x = c("lon", "lat"), by.y = c("lon", "lat"), all.x = TRUE)
    lhs<-krigingGrid[complete.cases(krigingGrid), c("lat", "lon", "Value")]
    krigingGrid<-krigingGrid[!complete.cases(krigingGrid), c("lat", "lon", "Value")]
    
    ##transform longitude
    krigingGrid$lon<-krigingGrid$lon - 360
    lhs$lon<-lhs$lon - 360
    ncdf4data$lon<-ncdf4data$lon - 360
    
    ## convert to geodata object
    geoKrigGrid<-as.geodata(krigingGrid)
    
    ##create empirical variogram
    evariog<-variog(geoKrigGrid, max.dist = 15, messages = FALSE)
  
    ## Fit five variogram models with optimized parameters
    fit1<-variofit(evariog, fix.nugget = FALSE, cov.model = "exponential",  messages = FALSE)
    fit2<-variofit(evariog, fix.nugget = FALSE, cov.model = "cauchy", messages = FALSE)
    fit3<-variofit(evariog, fix.nugget = FALSE, cov.model = "matern", messages = FALSE)
    fit4<-variofit(evariog, fix.nugget = FALSE, cov.model = "spherical", messages = FALSE)
    fit5<-variofit(evariog, fix.nugget = FALSE, cov.model = "cubic", messages = FALSE)
      
    ##  store them in a list
    fits<-list(fit1, fit2, fit3, fit4, fit5)
      
    ##  This loop is necessary becaus occaisionally the fitting process will produce a nugget such as -1e-15, which 
    ##  causes an error when trying to create a kriging control object. 
    for (i in seq_along(fits)){
      ## This is .1 becuase of precision problems in the logical operator, sometimes it doesn't realize it's less than zero
        if(fits[[i]]$nugget < .01) fits[[i]]$nugget<-0
    }
      
    ## experimental variogram fitting, add back to code when worked out. 
    ## system.time(smfit<-variofitsm(evariog, maxorder = 6, Rsqr = .98))
    
    ##create krige objects with each fit
    krigObj1<-krige.control(obj.model = fits[[1]])
    #krigObj2<-krige.control(obj.model = fits[[2]])
    krigObj3<-krige.control(obj.model = fits[[3]])
    krigObj4<-krige.control(obj.model = fits[[4]])
    krigObj5<-krige.control(obj.model = fits[[5]])
    
    ##store in list for selection of best model
    krigs<-list(krigObj1,  krigObj3, krigObj4, krigObj5)
    
    ##  perform Krigin with two different variogram models, one the best selected of five, the other always Cauchy (placeholder)
    ##  for standard model when code is ready
    krigPredict1<-krige.conv(geoKrigGrid, locations = lhs[,c("lat", "lon")], krige = krigObj1, output = list(messages = FALSE)) 
    #krigPredict2<-krige.conv(geoKrigGrid, locations = lhs[,c("lat", "lon")], krige = krigObj2, output = list(messages = FALSE))     
    krigPredict3<-krige.conv(geoKrigGrid, locations = lhs[,c("lat", "lon")], krige = krigObj3, output = list(messages = FALSE))  
    krigPredict4<-krige.conv(geoKrigGrid, locations = lhs[,c("lat", "lon")], krige = krigObj4, output = list(messages = FALSE))  
    krigPredict5<-krige.conv(geoKrigGrid, locations = lhs[,c("lat", "lon")], krige = krigObj5, output = list(messages = FALSE))  
       
    #record the prediction errors at the lhs points over oklahoma
    CV1<-sum((krigPredict1$predict - lhs$Value)^2)
    #CV2<-sum((krigPredict2$predict - lhs$Value)^2)
    CV3<-sum((krigPredict3$predict - lhs$Value)^2)
    CV4<-sum((krigPredict4$predict - lhs$Value)^2)
    CV5<-sum((krigPredict5$predict - lhs$Value)^2)
    
    ##get index of krig object with smallest cross validation SSRES
    CV<-c(CV1,  CV3, CV4, CV5)
    minCV<-which(CV == min(CV), arr.ind = TRUE)
    
    ##  perform kriging at the Mesonet locations using the model with lowest cross validation error
    MesoKrig<-krige.conv(geoKrigGrid, locations = stationlocs, krige = krigs[[minCV[1]]], output = list(messages = FALSE))  
    
    for ( i in 1:98) {
      
      ncvar_put(krigGEFfile, krigGEFSVAR, vals = MesoKrig$predict[i], start = c(i, d, t), count =c(1, 1, 1))
    
    }
   
    for ( i in 1:98) {
      
      ncvar_put(krigGEFfile, krigGEFSVARSE, vals = sqrt(MesoKrig$krige.var[i]), start = c(i, d, t), count =c(1, 1, 1))
      
    }

    
    Result<-data.frame(VarName, Day = t, SSRes = CV[minCV[1]])
    Result
  }##end d, t
  
  nc_close(krigGEFfile)
  
  Results
  
}

stopCluster(cl)




