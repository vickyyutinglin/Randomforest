# install.packages(c('readxl','caTools','randomForest','rgdal','ggplot2'))
# update.packages()
library(readxl)
library(caTools)
library(randomForest)
library(rgdal)
library(ggplot2)

setwd('C:\\Users\\USER.DESKTOP-ELH7F84\\OneDrive\\R project\\Randomforest')

# setwd('') # set your local path
stpi <- read_excel("Stylophora_pistillata.xlsx") # import Stylophora pistillata cover and environmental data from 209 transects
stpi <- stpi[ ,-c(1,5,6)] # remove the column with transect, min SST and max SST information
cord <- read_excel("cord.xlsx") # import the coordinate of sampling transects
taiwan <- readOGR('TWN_adm0.shp') # import the OGR file of Taiwan map 

knitr::opts_chunk$set(cache = T)

sample.stpi <- caTools::sample.split(stpi$SP, SplitRatio = .7)
stpi.train <- subset(stpi, sample.stpi == TRUE)
stpi.test <- subset(stpi, sample.stpi == FALSE)
dim(stpi.train)
dim(stpi.test)

set.seed(17) # set a random seed to make the model reproducible
rf.stpi <- randomForest(SP ~ .,data=stpi.train, importance = T) # build the random forest model
rf.stpi

set.seed(17) # set a random seed to make the model reproducible
plot(rf.stpi) # check the mean square error of models

which.min(rf.stpi$mse)

rf.stpi.min <- randomForest(SP ~ .,data=stpi.train, importance = T, ntree=21)
rf.stpi.min

rf.stpi.min.2 <- randomForest(SP ~ .,data=stpi.train, importance = T,ntree=21, mtry=2)
rf.stpi.min.2
rf.stpi.min.3 <- randomForest(SP ~ .,data=stpi.train, importance = T,ntree=21, mtry=3)
rf.stpi.min.3

varImpPlot(rf.stpi.min)
importance(rf.stpi.min)

pred.stpi <- predict(rf.stpi.min, newdata=stpi.test)
rf.resid.stpi <- pred.stpi - stpi.test$SP # the residuals
plot(stpi.test$SP,rf.resid.stpi, pch=18,ylab="residuals (predicted-observed)",   xlab="observed",col="blue3") # residual plot
abline(h=0,col="red4",lty=2)
mean(rf.resid.stpi^2) # mean squared error

par(mfrow=c(1,3)) 
partialPlot(rf.stpi.min, as.data.frame(stpi.train), light)
partialPlot(rf.stpi.min, as.data.frame(stpi.train), SST)
partialPlot(rf.stpi.min, as.data.frame(stpi.train), Chl_a)
dev.off()

cord <- cbind(cord,stpi$SP)
colnames(cord) <- c('Tr','Lat',"Lon",'Cover')
lat <- as.factor(cord$Lat)
lon <- as.factor(cord$Lon)
coords <- SpatialPoints(cord[, c("Lon", "Lat")]) 
stpi_cord <- SpatialPointsDataFrame(coords, cord) 
proj4string(stpi_cord) <- CRS("+proj=longlat +no_defs +ellps=WGS84", doCheckCRSArgs = F) 
head(stpi_cord)
class(stpi_cord)

taiwan_map <- fortify(taiwan)
Taiwan.map <- ggplot() +
  geom_polygon(data= taiwan_map, aes(x= long, y= lat, group= group), fill='light grey') +
  scale_x_continuous(limits = c(119, 123)) + # set x and y limt
  scale_y_continuous(limits = c(21.5, 25.5)) +
  theme_minimal()+  # background    
  coord_fixed(1.1)+ # fix the coordinate
  geom_point(cord, mapping = aes(x= Lon, y= Lat, size = Cover),alpha=0.2, color="blue") + # add the location of sampling transects
  theme(legend.position="right")

Taiwan.map+ scale_size_continuous(name = 'Cover (m^2)' , range = c(0, 15),breaks = c(0,5,10,15,20), labels = c('0.01','0.05','0.1','0.15','0.2')) # adjust the size of bubble plot
