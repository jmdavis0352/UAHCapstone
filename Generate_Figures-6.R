##  UAH MS OR Capstone research
##  Kaggle AMS Solar Enger Prediction Competition
##
##  Jonathan Davis 2014-03-15
##
##  Script 6: create images and plots for paper


library(maps)
library(scales)
library(ggplot2)
library(lhs)
library(geoR)
library(RColorBrewer)
library(grid)
library(rpart)
library(gbm)
library(rpart.plot)

##set working directory and load workspace
setwd("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/")
setwd("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/")


load("E:/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")
load("/media/jmdavis/StorageDrive/Dropbox/Kaggle/Solar Energy/R_ProjectSolarEnergy/Derived Data/PreparedData1.Rdata")


############    Figure 1 Gefs grid and mesonet locations   ###################
png("Derived Data/Figures/Fig1GEFSgrid.png",height = 600, width = 1066)
  par(mar = c(0.8, 0.5, 0.8, 0.5))
  map('county', 'oklahoma', ylim = c(31, 39), xlim = c(-106, -91) )
  points(expand.grid(Longs-360, Lats), pch = 19, col = alpha('black', .5), cex = 1.6 )
  points(stationinfo$elon, stationinfo$nlat, col = 'red', pch = 3, cex = .7)
legend("bottomleft", cex = 1.5, pch = c(3, 19), col = c('red', alpha('black', .5)), legend = c("Mesonet station", "GEFS grid"))

dev.off()


############    Figure 2 Solar irradiance  ###################
png("./Derived Data/Figures/Fig2ACME.png",height = 600, width = 1266)
par(mar = c(4.1,5,1,.1))
plot(train$ACME[1:2000], pch = 19, cex = 1.6, col = alpha("red", .4),cex.axis = 1.8, cex.lab = 2, xlab = "Day", ylab = expression(Jm^-2))
points(train$TAHL[1:2000], pch = 19, cex = 1.6, col = alpha("blue", .2))
dev.off()


############    Figure 3 Variogram and cross validation   ###################
load("./Derived Data/Figures/VariogramCompare.RData")
load("Derived Data/Figures/VariogramCompare.RData")

xy<-expand.grid(254:269, 31:39)
lhs<-optimumLHS(10, 2,eps = .05 ) ##generate design
lhs[,1]<-ceiling(lhs[,1]/(1/10)) + 256 ##convert to discrete longitude values
lhs[,2]<-ceiling(lhs[,2]/(1/4)) + 33 ##convert to discrete lattitude values
lhs[,1]<-lhs[,1] - 360
##plot the points to be removed
c<-2
png("./Derived Data/Figures/lhsExperiment.png", height = 700, width = 1400 )
  par(mar = c(4.5, 5.5, 2, 1), mfrow = c(1,2),xpd = NA)
  plot(evariog, cex = 1.4, cex.lab = c, cex.axis = c, xlab = "h", ylab = expression(gamma[e](h)))
  lines(fit1, col = 'black')## the exponential model
  lines(fit3, lty = 2, col = 'green') ## the linear model
  lines(fit4, lty = 3, col ='red') ## the spherical model
  lines(fit5, lty = 4, col = 'blue') ## cubic model
  text(0,max(evariog$v), "A) Variogram models", cex = c, adj = 0)
  legend("bottomright", pch = c(1,NA, NA,  NA, NA), cex =c, lty = c(NA,1,2,3,4), col = c('black', 'black', 'green', 'red', 'blue'),
         legend = c("Empirical variogram", "Exponential", "Matern", "Spherical", "Cubic") )


  map('county', 'oklahoma', ylim = c(31, 39), xlim = c(-106, -92) )
  points(expand.grid(Longs-360, Lats), pch = 19, col = alpha('black', .5), cex = 1.6 )
  points(lhs[,1], lhs[,2], pch = 19, col = 'white', cex = 1.8)
  points(lhs[,1], lhs[,2], pch = 10,cex = 1.6, col = 'red')
  text(-106, 39.9, "B) Data removed for variogram model cross validation", cex = c, adj = 0)
  legend("bottomright", cex = c,pch = c(19, 10, NA), col = c(alpha('black', .5), 'red', 'white'), legend = c("GEFS Grid", "Removed for cross validation", ""))
dev.off()



####  plot example of a regression tree #######
png("./Derived Data/Figures/RegTree.png", height = 1100, width = 1800, pointsize = 15)
rp<-rpart(form, MesonetC[MesonetC[,"Station"]=="ACME", ], method = "anova")
prp(rp)
dev.off()


################# Figure N Model Tuning  #################################
load("./Derived Data/Figures/GBMTuning3.RData")
tsize = 18
spectral<-brewer.pal(11,"RdYlBu")
png("Derived Data/Figures/GBMTuning.png", height = 1100, width = 1400, pointsize = 15)
g<-ggplot(gbm.Models, aes(x = Aggregated, y = Loss)) + geom_point(aes(shape = factor(shrinkage), colour = runtime), size = 4.5) + 
  scale_colour_gradientn(colours = rainbow(10), name = "Run Time (s)") + facet_grid(CVFolds~InteractionDepth, labeller = label_both) + 
  theme(axis.text.x = element_text(colour = 'grey30', size = tsize), axis.text.y = element_text(colour = "grey30", size = tsize), 
        axis.title.y = element_text(vjust = .01, size = tsize+2), strip.text.x = element_text(size = tsize), strip.text.y = element_text(size = tsize), 
        axis.title.x = element_text(size = tsize +2), legend.text = element_text(size = tsize), 
        legend.title = element_text(size = tsize), plot.margin=unit(c(1,3,2,6),"mm")) + 
  guides(shape = guide_legend(title = "Learning Rate", title.theme = element_text(size = tsize, face = 'bold',angle = 0))) +
  xlab("Number of Models Aggregated") + ylab("Mean Absolute Error") + 
  geom_point(data = gbm.Models[gbm.Models[, "shrinkage"]==.01,], aes(x= Aggregated, y = Loss), shape = 1, size = 4.5) +
  geom_point(data = gbm.Models[gbm.Models[, "shrinkage"]==.03,], aes(x= Aggregated, y = Loss), shape = 2, size = 4.5) +
  geom_point(data = gbm.Models[gbm.Models[, "shrinkage"]==.05,], aes(x= Aggregated, y = Loss), shape = 0, size = 4.5)

print(g )
##gg<- geom_point(aes(x = Aggregated, y = Loss, shape = factor(shrinkage)), size = 4) +scale_shape(solid = FALSE)


dev.off()

#########relative influence

#   load
rI<-summary(allModel2)

png("./Derived Data/Figures/RelInfluence.png", height = 600, width = 1400, pointsize = 15)
  par(las = 1, mar = c(5, 9, .1, 1), cex.lab = 1.6)
  barplot(rI$rel.inf[10:1], cex.names = 1.6, cex.axis = 1.4, horiz = TRUE, names.arg = rI$var[10:1], xlab = "Relative Influence")
dev.off()


