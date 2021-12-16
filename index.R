# install.packages(c('readxl','randomForest','rgdal','ggplot2','sf','tmap','sp'))
# update.packages()
library(readxl)
library(randomForest)
library(rgdal)
library(ggplot2)
library(sf)
library(tmap)
library(sp)

setwd('C:\\Users\\USER.DESKTOP-ELH7F84\\OneDrive\\R project\\Randomforest')

# setwd('') # set your local path
stpi.data <- read_excel("Stylophora_pistillata.xlsx") # import Stylophora pistillata cover and environmental data from 209 transects
stpi <- stpi.data[ ,-c(1,2)] # remove the column with transect information
cord <- read_excel("cord.xlsx") # import the coordinate of sampling transects
taiwan <- readOGR('TWN_adm0.shp') # import the OGR file of Taiwan map 
pred.env <- read.csv("env.csv") # import the coordinate of sampling transects
pred.env <- pred.env[,-1] 

knitr::opts_chunk$set(cache = T)

set.seed(1) # set a random seed to make the model reproducible
rf.stpi <- randomForest(SP ~ .,data=stpi, importance = T) # build the random forest model
rf.stpi

plot(rf.stpi) # check the mean square error of models

which.min(rf.stpi$mse)

rf.stpi.min <- randomForest(SP ~ ., data=stpi, importance = T, ntree=49)
rf.stpi.min

rf.stpi.min.2 <- randomForest(SP ~ ., data=stpi, importance = T, ntree=49, mtry=2)
rf.stpi.min.2
rf.stpi.min.3 <- randomForest(SP ~ ., data=stpi, importance = T, ntree=49, mtry=3)
rf.stpi.min.3

pred.stpi <- predict(rf.stpi.min, newdata=stpi)
rf.resid.stpi <- pred.stpi - stpi$SP # the residuals
plot(stpi$SP,rf.resid.stpi, pch=18,ylab="residuals (predicted-observed)",   xlab="observed",col="blue3") # residual plot
abline(h=0,col="red4",lty=2)

mean(rf.resid.stpi^2) # mean squared error

pred.env.5 <- data.frame(pred.env$light.5,pred.env$SST,pred.env$chl)
colnames(pred.env.5) <- c('light','Chl_a','SST')
pred.5 <- predict(rf.stpi.min, newdata=pred.env.5)

pred.env.10 <- data.frame(pred.env$light.10,pred.env$SST,pred.env$chl)
colnames(pred.env.10) <- c('light','Chl_a','SST')
pred.10 <- predict(rf.stpi.min, newdata=pred.env.10)

pred.env.15 <- data.frame(pred.env$light.15,pred.env$SST,pred.env$chl)
colnames(pred.env.15) <- c('light','Chl_a','SST')
pred.15 <- predict(rf.stpi.min, newdata=pred.env.15)

pred.env.40 <- data.frame(pred.env$light.40,pred.env$SST,pred.env$chl)
colnames(pred.env.40) <- c('light','Chl_a','SST')
pred.40 <- predict(rf.stpi.min, newdata=pred.env.40)

summary(pred.5)
summary(pred.10)
summary(pred.15)
summary(pred.40)

varImpPlot(rf.stpi.min)
importance(rf.stpi.min)

par(mfrow=c(1,3)) 
partialPlot(rf.stpi.min, as.data.frame(stpi), light)
partialPlot(rf.stpi.min, as.data.frame(stpi), SST)
partialPlot(rf.stpi.min, as.data.frame(stpi), Chl_a)
dev.off()

taiwan.sf <- st_as_sf(taiwan) # tranform the sp file to sf file 
cord.all <- cbind(cord, stpi.data, pred.stpi) # create a matrix with coordinate, observed, predicted cover for the sampling sites

# creat a function to match the sampling sites and the grid system I created
# I can specify the matix (with coordinate, depth, observed, predicted coral cover 
# information), the depths and the grid resolution (degree) I would like to use
grid.match.depth <- function(matrix, depth, resolution){
  cover <- subset(matrix , Depth==depth) # subset the matrix with the depth
  lat <- as.factor(cover$Lat) # setup the latitude as factor
  lon <- as.factor(cover$Lon) # setup the longitude as factor
  coords <- SpatialPoints(cover[, c("Lon", "Lat")]) # transfer the latitude & longitude into spatialpoint (for map)
  proj4string(coords) <- '+proj=longlat +datum=WGS84 +no_defs' # assign the coordinate projection system to the same one with the Taiwanese map
  
  # creat a grid system
  bb <- bbox(taiwan) # set up the boundary of this system to the same with the Taiwanese map
  cs <- c(1, 1)*resolution  # cell size (degree) : using the argument "resolution" to specify
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd) # cell topology
  sp_grd <- SpatialGridDataFrame(grd, data=data.frame(id=1:prod(cd)), proj4string='+proj=longlat +datum=WGS84 +no_defs') # create a sp data frame 
                                 
  grid.over <- cbind(over(coords, sp_grd),cover$SP ,cover$pred.stpi) # create a matrix with observed coral cover & predicted cover
  colnames(grid.over) <- c('id','obs.cover','pred.cover') # remane the matrix
  grid.over.mean <- aggregate(obs.cover~id, grid.over,mean) # at least 5 transects sampled in the same grid, average the observed cover of those transects
  loc <- as.numeric(grid.over.mean$id) # setup the id of the grids as numerical 
  id <- rep(1, prod(cd)) # setup a id vector with 1
  over <- replace(id,loc,grid.over.mean$obs.cover) # replace the ids with coral cover information with the mean observed coral cover
  grid.over.mean.pred <- aggregate(pred.cover~id, grid.over,mean) # at least 5 transects sampled in the same grid, average the predicted cover of those transects
  loc.pred <- as.numeric(grid.over.mean.pred$id) # setup the id of the grids as numerical 
  id.pred <- rep(1, prod(cd)) # setup a id vector with 1
  over.pred <- replace(id.pred,loc.pred,grid.over.mean.pred$pred.cover) # replace the ids with coral cover information with the mean predicted coral cover
  sp.grd <- SpatialGridDataFrame(grd,data=data.frame(observed.cover = over, predicted.cover = over.pred), # create a spatial data frame with the observed and predicted cover
                                 proj4string='+proj=longlat +datum=WGS84 +no_defs')
  
  return(list(grid.over,summary(sp.grd),sp.grd)) # return the result
}

layer.5m <- grid.match.depth(cord.all,5,0.1) # extract the information of 5m
sp.grd.5 <- layer.5m[[3]][1] # assign the 5m observed cover to the matrix
sp.grd.5.pred <- layer.5m[[3]][2] # assign the 5m predicted cover to the matrix
layer.10m <- grid.match.depth(cord.all,10,0.1)  # extract the information of 10m
sp.grd.10 <- layer.10m[[3]][1] # assign the 10m observed cover to the matrix
sp.grd.10.pred <- layer.10m[[3]][2] # assign the 10m predicted cover to the matrix
layer.15m <- grid.match.depth(cord.all,15,0.1)  # extract the information of 15m
sp.grd.15 <- layer.15m[[3]][1] # assign the 15m observed cover to the matrix
sp.grd.15.pred <- layer.15m[[3]][2] # assign the 15m predicted cover to the matrix
layer.40m <- grid.match.depth(cord.all,40,0.1)  # extract the information of 40m
sp.grd.40 <- layer.40m[[3]][1] # assign the 40m observed cover to the matrix
sp.grd.40.pred <- layer.40m[[3]][2] # assign the 40m predicted cover to the matrix

bb <- bbox(taiwan) # set up the boundary of this system to the same with the Taiwanese map
cs <- c(1, 1)*0.1  # cell size 0.1 x 0.1 degree 
cc <- bb[, 1] + (cs/2)  # cell offset
cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd) # cell topology
sp_grd <- SpatialGridDataFrame(grd,data=data.frame(id=1:prod(cd)), proj4string='+proj=longlat +datum=WGS84 +no_defs')

sp.grd.all <- SpatialGridDataFrame(grd,data=data.frame(coral.cover.5m = sp.grd.5$observed.cover, # create a spatial data frame with the observed cover for mapping
                                                       coral.cover.10m = sp.grd.10$observed.cover,
                                                       coral.cover.15m = sp.grd.15$observed.cover,
                                                       coral.cover.40m = sp.grd.40$observed.cover),
                                   proj4string='+proj=longlat +datum=WGS84 +no_defs')
sp.grd.all.sf <- st_as_sf(sp.grd.all) # tranform the sp file to the sf file

sp.grd.pred.all <- SpatialGridDataFrame(grd,data=data.frame(coral.cover.5m = sp.grd.5.pred$predicted.cover, # create a spatial data frame with the predicted cover in sampling sites for mapping
                                                            coral.cover.10m = sp.grd.10.pred$predicted.cover,
                                                            coral.cover.15m = sp.grd.15.pred$predicted.cover,
                                                            coral.cover.40m = sp.grd.40.pred$predicted.cover),
                                        proj4string='+proj=longlat +datum=WGS84 +no_defs')
sp.grd.all.pred.sf <- st_as_sf(sp.grd.pred.all) # tranform the sp file to the sf file
# draw the maps
col.pan <- c("lightblue","blue","green", "yellow","orange","red","white") # setup the color panel for mapping the cover
tmap_mode("plot")

tm_shape(sp.grd.all.sf) +
  tm_symbols(c("coral.cover.5m","coral.cover.10m","coral.cover.15m","coral.cover.40m"),
             size=1,breaks = c(seq(0,0.12,0.02),1),shape=15, 
             palette = col.pan, title.size = 2, border.col =NA,
             labels = c("0-0.02", "0.02-0.04", "0.04-0.06","0.06-0.08","0.08-0.10","0.10-0.12","No cover data"))+
tm_shape(taiwan.sf) +
  tm_polygons("Islands",legend.show = F, alpha =.3,)
    
tm_shape(sp.grd.all.pred.sf) +
  tm_symbols(c("coral.cover.5m","coral.cover.10m","coral.cover.15m","coral.cover.40m"),
             size=1,breaks = c(seq(0,0.12,0.02),1),shape=15, 
             palette = col.pan,title.size = 2, border.col =NA,
             labels = c("0-0.02", "0.02-0.04", "0.04-0.06","0.06-0.08","0.08-0.10","0.10-0.12","No cover data"))+
tm_shape(taiwan.sf) +
  tm_polygons("Islands",legend.show = F, alpha =.3,)

# grid with coral cover prediction
pred.5.map <- pred.5 ; pred.10.map <- pred.10
pred.15.map <- pred.15 ; pred.40.map <- pred.40
pred.5.map[is.na(pred.5.map)] <- 1 # replace the na with 1
pred.10.map[is.na(pred.10.map)] <- 1
pred.15.map[is.na(pred.15.map)] <- 1
pred.40.map[is.na(pred.40.map)] <- 1
sp.grd.pred.without <- SpatialGridDataFrame(grd,data=data.frame(Predicted.coral.cover.at.5m = pred.5.map, # create a spatial data frame with the predicted cover of all grids for mapping
                                                                Predicted.coral.cover.at.10m = pred.10.map,
                                                                Predicted.coral.cover.at.15m = pred.15.map,
                                                                Predicted.coral.cover.at.40m = pred.40.map),
                                            proj4string='+proj=longlat +datum=WGS84 +no_defs')
sp.grd.pred.without.sf <- st_as_sf(sp.grd.pred.without) # tranform the sp file to the sf file

# draw the maps for the predicted coral cover around Taiwan
tm_shape(sp.grd.pred.without.sf) +
  tm_symbols(c("Predicted.coral.cover.at.5m","Predicted.coral.cover.at.10m", "Predicted.coral.cover.at.15m","Predicted.coral.cover.at.40m"),
             size=1,breaks = c(seq(0,0.12,0.02),1),
             palette = col.pan, title.size = 5,border.col =NA,shape=15, 
             labels = c("0-0.02", "0.02-0.04", "0.04-0.06","0.06-0.8","0.08-0.10","0.10-0.12","Missing data")) +
tm_shape(taiwan.sf) +
  tm_polygons("Islands",legend.show = F,alpha =.3)
