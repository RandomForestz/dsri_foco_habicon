#Running SS Optim for pallida on Pegasus

library(ResistanceGA)

#set up Julia
JULIA_HOME <- "/nethome/ccm151/src/julia-1.5.2/bin" 
JuliaCall::julia_setup(JULIA_HOME)

setwd("/scratch/projects/amphib_demo/caitlin/pond_turtle")

#load in points
pall <- read_csv("data/pall_us_pop.csv")

#make spatial and change projection

#convert to spatial points and filter out unique cords
pall_sp <- SpatialPointsDataFrame(pall[,c('longitude', 'latitude'),], data = pall)

proj4string(pall_sp) <- CRS("+proj=longlat +datum=WGS84")

newprj <-  crs("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs")

pall_sp <- spTransform(pall_sp, newprj) #need to get in same projection

# filter out unique sites
pall_pop <-  pall_sp[which(!duplicated(pall_sp$site_number)), ] %>% 
  SpatialPoints(.)

#load in genetic distance matrix
pall_fst <- read.csv("data/pall_fst_us.csv")

row.names(pall_fst) <- NULL

## GA.prep and jl.prep

#make directory for outputs
dir.create(file.path("data/", "pall_SSOptim_results"))
#name that directory
results_dir <- "data/pall_SSOptim_results/"
#name directory where surfaces are
rast_dir <- "data/pall_surfaces/"


GA.inputs <- GA.prep(ASCII.dir = rast_dir, 
                     Results.dir = results_dir, parallel = TRUE)

jl.inputs <- jl.prep(n.Pops = length(pall_pop), 
                     response = lower(pall_fst), 
                     CS_Point.File = pall_pop, JULIA_HOME = JULIA_HOME,
                     parallel = TRUE, silent = FALSE)


## Single Surface Optimization

jl.optim <- SS_optim(jl.inputs = jl.inputs, GA.inputs = GA.inputs)

