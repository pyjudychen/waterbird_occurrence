# Data preparation----
## Load packages----
library(rgbif)
library(maptools)
# install.packages('ENMeval')
library(ENMeval)
library(raster)
library(dismo)
# install.packages('rJava')

Sys.setenv('JAVA_HOME'="C:/Program Files/Java/jdk-19/")
library(rJava)


fields <- c('name', 'basisOfRecord', 'scientificName', 'county','state', 'catalogNumber', 'institutionCode', 
            'decimalLatitude', 'decimalLongitude', 'year', 'month', 
            'stateProvince', 'country', 'habitat', 'family', 'genus', 'locality')

spp <- occ_search(scientificName = c("Anas acuta"), fields = fields, institutionCode = 'CLO', stateProvince = 'California', 
                  limit = 20000, year = 2021)

View(spp)
View(spp$data)
spp.data <- spp$data
str(spp.data)

spp.data$lat <- spp.data$decimalLatitude
spp.data$lon <- spp.data$decimalLongitude

aacuta <- spp.data

data(wrld_simpl)
plot(wrld_simpl, xlim=c(-130,-80), ylim=c(20, 50), axes=TRUE, col="light yellow")+
  box()+
  points(aacuta$lon, aacuta$lat, col = 'red', cex = 0.75)


##Data cleaning----
### Data without coordinates---- 
lonzero = subset(aacuta, lon == 0) 

### Duplicate records----
dups <- duplicated(aacuta[, c('lon', 'lat')])
sum(dups)

aacuta <- aacuta[!dups, ]
head(aacuta)
str(aacuta)

print(cat('Records remain: ', as.character(nrow(aacuta))))

## Save the cleaned file----
write.csv(aacuta, "C:\\Users\\judy8\\Box\\PhD_UCD\\Data\\aacuta_clean.csv")
aacuta <- read.csv('C:\\Users\\judy8\\Box\\PhD_UCD\\Data\\aacuta_clean.csv')
str(aacuta)

# Environmental data----
envs <- raster::getData(name = "worldclim", var = "bio", res = 10, lat = c(20, 50), lon = c(-130,-80))
envRes <- 10
if (envRes == 0.5) {
  i <- grep('_', names(envs))
  editNames <- sapply(strsplit(names(envs)[i], '_'), function(x) x[1])
  names(envs)[i] <- editNames
}
class(envs)

i <- grep('bio[0-9]$', names(envs))
editNames <- paste('bio', sapply(strsplit(names(envs)[i], 'bio'), function(x) x[2]), sep='0')
names(envs)[i] <- editNames

envs <- envs[[c('bio01', 'bio04', 'bio05', 'bio06', 'bio12', 'bio16', 'bio17')]]
plot(envs, xlim = c(-130,-80), ylim = c(20, 50))

# BIO1 = Annual Mean Temperature
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month**
# BIO12 = Annual Precipitation
# BIO16 = Precipitation of Wettest Quarter**
# BIO17 = Precipitation of Driest Quarter**

data(wrld_simpl)
str(aacuta)

plot(envs, 1, xlim = c(-130,-80), ylim = c(20, 50)) 
plot(wrld_simpl, add = TRUE)+
  # mtext('Annual Mean Temperature')
  points(x = aacuta$lon, y = aacuta$lat, col = 'blue')

colnames(aacuta)  

## Extract values from rasters----
# extract environmental values at aacuta grid cells
locs.vals <- raster::extract(envs[[1]], aacuta[, c('lon', 'lat')])
# extract environmental values at occ grid cells
aacuta <- aacuta[!is.na(locs.vals), ]  
print(cat("Records remain: ", as.character(nrow(aacuta))))

xmin <- min(aacuta$lon)
xmax <- max(aacuta$lon)
ymin <- min(aacuta$lat)
ymax <- max(aacuta$lat)
bb <- matrix(c(xmin, xmin, xmax, xmax, xmin, ymin, ymax, ymax, ymin, ymin), ncol=2)
bgExt <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(bb)), 1)))
bgExt <- rgeos::gBuffer(bgExt, width = 0.5)

# crop the environmental rasters by the background extent shape
envsBgCrop <- raster::crop(envs, bgExt)
# mask the background extent shape from the cropped raster
envsBgMsk <- raster::mask(envsBgCrop, bgExt)
# sample random background points
bg.xy <- dismo::randomPoints(envsBgMsk, 5000)
# convert matrix output to data frame
bg.xy <- as.data.frame(bg.xy)
print(cat("Number of background points: ", as.character(nrow(bg.xy))))

#################################################################
aacuta <- aacuta[, c(17,16)]
presvals <- extract(envs, aacuta)
set.seed(0)
backgr <- randomPoints(envs, 500)
absvals <- extract(envs, backgr)
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))

sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

head(sdmdata)
tail(sdmdata)

pairs(sdmdata[, 2:5], cex = 0.1)

saveRDS(sdmdata, 'C:\\Users\\judy8\\Box\\PhD_UCD\\Results\\sdm.Rds')
saveRDS(presvals, "C:\\Users\\judy8\\Box\\PhD_UCD\\Results\\pvals.Rds")
#################################################################

# Prepare Project Niche Model Enveronment Background data---- 
projCoords <- data.frame(x = c(-80, -80, -130, -130, -80), y = c(60, 20, 20, 60, 60))
projPoly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(projCoords)), ID=1)))
plot(projPoly)
plot(wrld_simpl, add = TRUE)

# crop the environmental rasters by the background extent shape
predsProj <- raster::crop(envs, projPoly)
# mask the background extent shape from the cropped raster
predsProj <- raster::mask(predsProj, projPoly)

#############################################################################
# Model fitting, prediction, and evaluation----
## Model fitting----
sdmdata <- readRDS('C:\\Users\\judy8\\Box\\PhD_UCD\\Results\\sdm.Rds')
presvals <- readRDS('C:\\Users\\judy8\\Box\\PhD_UCD\\Results\\pvals.Rds')

m1 <- glm(pb ~ bio01 + bio05 + bio12, data = sdmdata)
summary(m1)

m2 <- glm(pb ~., data = sdmdata)
summary(m2)

m3 <- glm(pb ~ bio05 + bio06 + bio12 + bio16, data = sdmdata)
summary(m3)

# most models only take presence points
bc <- bioclim(presvals[, c('bio01', 'bio05', 'bio12')])
class(bc)

bc
pairs(bc)


## Model prediction----
bio01 = c(40, 150, 200)
bio05 = c(60, 115, 290)
bio12 = c(600, 1600, 1700)
pd = data.frame(cbind(bio01, bio05, bio12))
pd

predict(m1, pd)
predict(bc, pd)
response(bc)

names(envs)
p <- predict(envs, m2)
plot(p)

## Model evaluation----
p <- rnorm(50, mean = 0.7, sd = 0.3)
a <- rnorm(50, mean = 0.4, sd = 0.4)

par(mfrow = c(1, 2))
plot(sort(p), col = 'tomato', pch = 21)+
  points(sort(a), col = 'navy', pch = 24)+
  legend(1, 0.95*max(a, p), c('presence', 'absence'), 
         pch = c(21, 24), col = c('tomato', 'navy'))
comb <- c(p, a)
group <- c(rep('presence', length(p)), rep('absence', length(a)))
boxplot(comb ~ group, col = c('navy', 'tomato'))

###############################################################################

# Build and Evaluate Niche Model----
## dismo::maxent----
occs.xy <- aacuta[c('lon', 'lat')]
group.data <- ENMeval::get.randomkfold(occ=occs.xy, bg=bg.xy, kfolds=3)

# pull out the occurrence and background partition group numbers from the list
occs.grp <- group.data[[1]]
bg.grp <- group.data[[2]]

fold <- kfold(occs.xy, k=5)
occtest <- occs.xy[fold == 1, ]
occtrain <- occs.xy[fold != 1, ]

# Run Maxent Model & see the maxent results in a browser----
me <- maxent(envsBgMsk, data.frame(occtrain), args=c("-J", "-P"))
me
response(me)
plot(me)

e1 <- evaluate(me, p = occtest, a = bg.xy, x = envsBgMsk)
plot(e1, 'ROC')+
  mtext(substitute(paste('Sensitivity vs. 1-specificity for ', italic('Anas acuta'))))
  
threshold(e1)

r <- predict(me, envsBgMsk)
plot(r)+
  points(occs.xy, col = col.alpha('black', 0.2), pch = 16, cex = 0.5)+
  mtext(substitute(paste('Suitability of locations in 2021 for ', italic('Anas acuta'))))