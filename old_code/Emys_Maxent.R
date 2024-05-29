# Emys Maxent model

library(rgdal)
library(rJava)
options(java.parameters = "-Xmx14g" )
library(dismo)

rasterOptions(tmpdir = 'temp/', progress = "text", maxmemory = 14e9)
removeTmpFiles(h = 17)


##MARMORATA --------------------------------------------------------------


#List of species (just one in this case but makes reusing code easier)
s <- c("Emys marmorata")
specs <- list(s)

#presence data
###Check class to make sure spatial points data frame

marm <- read.csv("data/pondTurt/marm_pop.csv")

#make spatial and change projection

#convert to spatial points and filter out unique cords
coordinates(marm) <- c("longitude", "latitude")

proj4string(marm) <- CRS("+proj=longlat +datum=WGS84")

newprj <-  crs("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs")

marm_sp <- spTransform(marm, newprj) #need to get in same projection

# filter out unique sites
marm_pop <-  marm_sp[which(!duplicated(marm_sp$site_number)), ]


#background points

bg_marm <- read.csv("data/pondTurt/bg_marm.csv")
coordinates(bg_marm) <- c("decimalLongitude", "decimalLatitude")

proj4string(bg_marm) <- CRS("+proj=longlat +datum=WGS84")
bg_marm <- spTransform(bg_marm, newprj)



#Environmental predictors
marm_preds <- brick("D:/PondTurtle/marm_preds.tif")

#assign proper names to predictors
names(marm_preds) <- c("wetlands", "maxtemp", "mintemp", "precip", "solar", 
                       "roads", "canopy", "developed", "elevation")


###########################################Current Model with all betas################################
#Set wd to output folder 
#dir.create("D:/PondTurtle/marm_maxent")

outpath <- "D:/PondTurtle/marm_maxent"
setwd(outpath)

#make sequence of betas to run loop for different beta values
betas <- array(dim = c(24, 2))
betas[,1] <- "betamultiplier="
betas[,2] <- c(seq(0.2,1, 0.2), seq(2, 20, 1))
betas <- data.frame(betas)
betas[,3] <- paste(betas[,1], betas[,2],sep = '')
names(betas)[3] <- 'setbeta'
betas[,2] <- as.character(betas[,2])

for (i in 1:length(specs)){    
  
  occ <- coordinates(marm_pop)
  
  bg <- coordinates(bg_marm)
  
  
  allbetas <- data.frame(array(dim = c(24, 7)))
  names(allbetas) <- c('species','beta','no.params','no.pres','loglik','AIC','AICc')
  allbetas[,1] <- specs[i]
  
  for (k in 1:nrow(betas)){
    
    
    xm <- maxent(x = marm_preds, p = occ, a = bg, args = c(betas[k,"setbeta"],"responsecurves=TRUE", "writebackgroundpredictions=TRUE"), path = outpath)
    
    save(xm, file = paste(specs[i], betas[k,"X2"], "xm.RData", sep = '_'))
    
    #get lambdas info out of xm
    lambdas <- data.frame(strsplit(xm@lambdas, ','))
    lambdas <- t(lambdas)
    lambdas <- data.frame(lambdas)
    rownames(lambdas) <- 1:nrow(lambdas)
    #vitalconstants <- lambdas[(nrow(lambdas)-3):nrow(lambdas),1:2]
    lambdas <- lambdas[-((nrow(lambdas)-3):nrow(lambdas)),]
    colnames(lambdas) <- c('variable','lambda_estimate','min','max')
    lambdas$lambda_estimate <- as.numeric(as.character(lambdas$lambda_estimate))
    
    #read and write predictions at presences
    pres.raw <- read.csv('species_samplePredictions.csv')
    write.csv(pres.raw,paste(specs[i], betas[k,"X2"], "pres.predictions.csv", sep = '_') )
    
    #calculate AICc
    no.params <- length(which(lambdas$lambda_estimate != 0))
    no.pres <- xm@results[1]
    loglik <- sum(log(pres.raw$Raw.prediction))
    AIC <- (2*no.params) -(2*loglik)
    AICc <- AIC + (((2*no.params)*(no.params+1))/(no.pres-no.params-1))
    
    
    allbetas[k,2:7] <-c(betas[k,"X2"], no.params, no.pres, loglik, AIC, AICc)
    
    
    
    print(k)
  }
  
  write.csv(allbetas, file = paste(specs[i],"allbetas.csv", sep="_"))
  print(i)
}


#make output data frame for beta results
bestbetas <- data.frame(array(dim = c(1, 5)))
names(bestbetas) <- c('species','no.pres','no.params default','no.params best','bestbeta')


for (i in 1:length(specs)){
  #load beta file
  inpath <- outpath
  f <- list.files(path = inpath, pattern = 'allbetas')
  setwd(inpath)
  data <- read.csv(f)
  
  
  bestbetas[i,1] <- as.character(data$species[1])
  bestbetas[i,2] <- data$no.pres[1]
  bestbetas[i,3] <- data[which(data$beta == 1),'no.params']
  
  if(max(data$no.pres) > min(data$no.params))
    bestbetas[i,5] <- min(data[which(data$AICc == min(data[which(data$no.params < data$no.pres&data$no.params >0),'AICc'])),'beta'])
  
  if((max(data$no.pres) > min(data$no.params)))
    bestbetas[i,4] <- data[which(data$beta == bestbetas[i,5]),'no.params']
  print(i)
}

save(bestbetas, file = 'bestbetas.RData')
write.csv(bestbetas, file = 'bestbetas.csv')

##################################Curent Model with  best beta######################

#set new directory for current best beta
dir.create("D:/PondTurtle/marm_maxent/best_beta")
outpath <- "D:/PondTurtle/marm_maxent/best_beta"
setwd(outpath)

for (i in 1:length(specs)){
  
  
  xm <- maxent(x = marm_preds, p = occ, a = bg,
               args = c(betas[which(betas$X2 == bestbetas[i,"bestbeta"]),"setbeta"],"responsecurves=TRUE", "writebackgroundpredictions=TRUE", "replicates=10"), 
               path = outpath)
  
  
  save(xm, file = paste(specs[i], bestbetas[i,"bestbeta"], "xm.RData", sep = '_'))
  
  write.csv(data.frame(xm@results[7:44,]), file = paste(specs[i],bestbetas[i,"bestbeta"],'variable_scores.csv', sep = '_'))
  
  #predict onto current rasters
  
  
  predict(xm, marm_preds, filename = paste(specs[i],'current_prediction.tif', sep = '_'), overwrite = TRUE, progress = 'text')
  
  print(i)
}

###CALCULATE AVERAGE OF PREDICT MAPS TO GET ONE SUITABILITY MAP
for (i in 1:length(specs)){
  
  
  f <- list.files(path = getwd(), pattern = 'current_prediction')
  ras <- lapply(f, raster)
  STACK1 <- stack(ras) 
  mean <- stackApply(STACK1, indices = rep(1,nlayers(STACK1)), fun = "mean", na.rm = T)
  writeRaster( x = mean, filename = paste(specs[i], 'current_suitability_map.tif', sep = '_'), overwrite = TRUE)
  
  
  print(i)
  
}


# PALLIDA ----------------------------------------------------------------


#List of species (just one in this case but makes reusing code easier)
s <- c("Emys pallida")
specs <- list(s)

#presence data
###Check class to make sure spatial points data frame

pall <- read.csv("data/pondTurt/pall_us_pop.csv")

#make spatial and change projection

#convert to spatial points and filter out unique cords
coordinates(pall) <- c("longitude", "latitude")

proj4string(pall) <- CRS("+proj=longlat +datum=WGS84")

newprj <-  crs("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs")

pall_sp <- spTransform(pall, newprj) #need to get in same projection

# filter out unique sites
pall_pop <-  pall_sp[which(!duplicated(pall_sp$site_number)), ]


#background points

bg_pall <- read.csv("data/pondTurt/bg_pall.csv")
coordinates(bg_pall) <- c("decimalLongitude", "decimalLatitude")

proj4string(bg_pall) <- CRS("+proj=longlat +datum=WGS84")
bg_pall <- spTransform(bg_pall, newprj)



#Environmental predictors
pall_preds <- brick("D:/PondTurtle/pall_preds.tif")

#assign proper names to predictors
names(pall_preds) <- c("wetlands", "maxtemp", "mintemp", "precip", "solar", 
                       "roads", "canopy", "developed", "elevation")


###########################################Current Model with all betas################################
#Set wd to output folder 
dir.create("D:/PondTurtle/pall_maxent")

outpath <- "D:/PondTurtle/pall_maxent"
setwd(outpath)

#make sequence of betas to run loop for different beta values
betas <- array(dim = c(24, 2))
betas[,1] <- "betamultiplier="
betas[,2] <- c(seq(0.2,1, 0.2), seq(2, 20, 1))
betas <- data.frame(betas)
betas[,3] <- paste(betas[,1], betas[,2],sep = '')
names(betas)[3] <- 'setbeta'
betas[,2] <- as.character(betas[,2])

for (i in 1:length(specs)){    
  
  occ <- coordinates(pall_pop)
  
  bg <- coordinates(bg_pall)
  
  
  allbetas <- data.frame(array(dim = c(24, 7)))
  names(allbetas) <- c('species','beta','no.params','no.pres','loglik','AIC','AICc')
  allbetas[,1] <- specs[i]
  
  for (k in 1:nrow(betas)){
    
    
    xm <- maxent(x = pall_preds, p = occ, a = bg, args = c(betas[k,"setbeta"],"responsecurves=TRUE", "writebackgroundpredictions=TRUE"), path = outpath)
    
    save(xm, file = paste(specs[i], betas[k,"X2"], "xm.RData", sep = '_'))
    
    #get lambdas info out of xm
    lambdas <- data.frame(strsplit(xm@lambdas, ','))
    lambdas <- t(lambdas)
    lambdas <- data.frame(lambdas)
    rownames(lambdas) <- 1:nrow(lambdas)
    #vitalconstants <- lambdas[(nrow(lambdas)-3):nrow(lambdas),1:2]
    lambdas <- lambdas[-((nrow(lambdas)-3):nrow(lambdas)),]
    colnames(lambdas) <- c('variable','lambda_estimate','min','max')
    lambdas$lambda_estimate <- as.numeric(as.character(lambdas$lambda_estimate))
    
    #read and write predictions at presences
    pres.raw <- read.csv('species_samplePredictions.csv')
    write.csv(pres.raw,paste(specs[i], betas[k,"X2"], "pres.predictions.csv", sep = '_') )
    
    #calculate AICc
    no.params <- length(which(lambdas$lambda_estimate != 0))
    no.pres <- xm@results[1]
    loglik <- sum(log(pres.raw$Raw.prediction))
    AIC <- (2*no.params) -(2*loglik)
    AICc <- AIC + (((2*no.params)*(no.params+1))/(no.pres-no.params-1))
    
    
    allbetas[k,2:7] <-c(betas[k,"X2"], no.params, no.pres, loglik, AIC, AICc)
    
    
    
    print(k)
  }
  
  write.csv(allbetas, file = paste(specs[i],"allbetas.csv", sep="_"))
  print(i)
}


#make output data frame for beta results
bestbetas <- data.frame(array(dim = c(1, 5)))
names(bestbetas) <- c('species','no.pres','no.params default','no.params best','bestbeta')


for (i in 1:length(specs)){
  #load beta file
  inpath <- outpath
  f <- list.files(path = inpath, pattern = 'allbetas')
  setwd(inpath)
  data <- read.csv(f)
  
  
  bestbetas[i,1] <- as.character(data$species[1])
  bestbetas[i,2] <- data$no.pres[1]
  bestbetas[i,3] <- data[which(data$beta == 1),'no.params']
  
  if(max(data$no.pres) > min(data$no.params))
    bestbetas[i,5] <- min(data[which(data$AICc == min(data[which(data$no.params < data$no.pres&data$no.params >0),'AICc'])),'beta'])
  
  if((max(data$no.pres) > min(data$no.params)))
    bestbetas[i,4] <- data[which(data$beta == bestbetas[i,5]),'no.params']
  print(i)
}

save(bestbetas, file = 'bestbetas.RData')
write.csv(bestbetas, file = 'bestbetas.csv')

##################################Curent Model with  best beta######################

#set new directory for current best beta
dir.create("D:/PondTurtle/pall_maxent/best_beta")
outpath <- "D:/PondTurtle/pall_maxent/best_beta"
setwd(outpath)

for (i in 1:length(specs)){
  
  
  xm <- maxent(x = pall_preds, p = occ, a = bg,
               args = c(betas[which(betas$X2 == bestbetas[i,"bestbeta"]),"setbeta"],"responsecurves=TRUE", "writebackgroundpredictions=TRUE", "replicates=10"), 
               path = outpath)
  
  
  save(xm, file = paste(specs[i], bestbetas[i,"bestbeta"], "xm.RData", sep = '_'))
  
  write.csv(data.frame(xm@results[7:44,]), file = paste(specs[i],bestbetas[i,"bestbeta"],'variable_scores.csv', sep = '_'))
  
  #predict onto current rasters
  
  
  predict(xm, pall_preds, filename = paste(specs[i],'current_prediction.tif', sep = '_'), overwrite = TRUE, progress = 'text')
  
  print(i)
}

###CALCULATE AVERAGE OF PREDICT MAPS TO GET ONE SUITABILITY MAP
for (i in 1:length(specs)){
  
  
  f <- list.files(path = getwd(), pattern = 'current_prediction')
  ras <- lapply(f, raster)
  STACK1 <- stack(ras) 
  mean <- stackApply(STACK1, indices = rep(1,nlayers(STACK1)), fun = "mean", na.rm = T)
  writeRaster( x = mean, filename = paste(specs[i], 'current_suitability_map.tif', sep = '_'), overwrite = TRUE)
  
  
  print(i)
  
}

#create inverse suitability maps to use as resistance for circuitscape
setwd("../../")

marm_suit <- raster("marm_maxent/best_beta/Emys marmorata_current_suitability_map.tif")

pall_suit <- raster("pall_maxent/best_beta/Emys pallida_current_suitability_map.tif")

marm_resist <- 1 - marm_suit
writeRaster(marm_resist, "marm_maxent/best_beta/Emys marmorata_resistance.tif")
writeRaster(marm_resist, "data/pondTurt/marm_resist.asc")
  
pall_resist <- 1 - pall_suit
writeRaster(pall_resist, "pall_maxent/best_beta/Emys pallida_resistance.tif")
writeRaster(pall_resist, filename = "~/Desktop/ConnectivityPrioritization/data/pondTurt/pall_resist.asc", overwrite = TRUE)
#Circuitscape does not use 0 resistance values, change to very small #
pall_resist[pall_resist == 0] <- 6e-08 #looked at next smallest value, saved this file

# CIRCUITSCAPE --------------------------------------------------------------------------

library(ResistanceGA)

JULIA_HOME <- "C:/Users/Caitlin/AppData/Local/Programs/Julia 1.5.3/bin/" 
JuliaCall::julia_setup(JULIA_HOME = JULIA_HOME) 

pall_pop <- SpatialPoints(pall_pop)
#rasterize points and save as ascii
pall_ascii <- rasterize(pall_pop, pall_resist)
writeRaster(pall_ascii, "~/Desktop/ConnectivityPrioritization/data/pondTurt/pall_points.asc")

marm_ascii <- SpatialPoints(marm_pop) %>% rasterize(., marm_resist)
writeRaster(marm_ascii, "data/pondTurt/marm_points.asc")

marm_CS <- Run_CS.jl(r = "~/Desktop/ConnectivityPrioritization/data/pondTurt/pall_resist.asc", CurrentMap = TRUE, output = "raster", EXPORT.dir = "~/Desktop/ConnectivityPrioritization/data/pondTurt/",
                     CS_Point.File = pall_pop, JULIA_HOME = JULIA_HOME, cholmod = FALSE)
#this wouldn't work, ran in Julia instead

#read in current map from circuitscape outputs
pall_current <- raster("data/pondTurt/pall_CS_cum_curmap.asc")

#crop pallida resitance to LA area and run circuitscape again


pall_resist <- raster("D:/PondTurtle/pall_maxent/best_beta/Emys pallida_resistance.tif")
plot(pall_resist)
points(pall_pop, pch = 16)

ext <- drawExtent()

pall_resist_LA <- crop(pall_resist, ext)

writeRaster(pall_resist_LA, "data/pondTurt/pall_resist_LA.tif", overwrite = TRUE)
writeRaster(pall_resist_LA, "data/pondTurt/pall_resist_LA.asc", overwrite = TRUE)

# clip points to LA regions

pall_pop_LA <- crop(pall_pop, extent(pall_resist_LA))

#save as csv and acii for circuitscape
pall_pop_LA <- SpatialPoints(pall_pop_LA)
#rasterize points and save as ascii
pall_ascii_LA <- rasterize(pall_pop_LA, pall_resist_LA)
writeRaster(pall_ascii_LA, "data/pondTurt/pall_points_LA.asc")
write.csv(pall_pop_LA, "data/pondTurt/pall_pop_LA.csv")
