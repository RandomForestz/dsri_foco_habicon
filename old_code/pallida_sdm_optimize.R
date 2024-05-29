# pallida SDM resistanceGA
## Set up julia

JULIA_HOME <- "C:/Users/Caitlin/AppData/Local/Programs/Julia 1.5.3/bin/" 
JuliaCall::julia_setup(JULIA_HOME = JULIA_HOME) #don't think need rebuild = true since run already


library(ResistanceGA)
library(rgdal)


## load in coordinate file
pall <- read.csv("data/pondTurt/pall_us_pop.csv")

#convert to spatial points and filter out unique cords
pall_sp <- SpatialPointsDataFrame(pall[,c('longitude', 'latitude'),], data = pall)

proj4string(pall_sp) <- CRS("+proj=longlat +datum=WGS84")

newprj <-  crs("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs")

pall_sp <- spTransform(pall_sp, newprj) #need to get in same projection

# filter out unique sites
pall_pop <-  pall_sp[which(!duplicated(pall_sp$site_number)), ] %>% 
  SpatialPoints(.)


## genetic dist file

pall_fst <- read.csv("data/pondTurt/pall_fst_us.csv")
pall_fst <- pall_fst[,-1]

#use chord dist
pall_chord <- read.csv("data/pondTurt/pall_chord.csv")
pall_chord <- pall_chord[,-1]


# read in sdm to be optimized
pall_sdm <- raster("D:/PondTurtle/pall_maxent/best_beta/Emys pallida_current_suitability_map.tif")

writeRaster(pall_sdm, "data/pondTurt/pall_sdm/pall_sdm.asc")


#make directory for outputs
dir.create("data/pondTurt/pallida_sdm_optim")
results_dir <- "data/pondTurt/pallida_sdm_optim/"
dir.create("data/pondTurt/pallida_sdm_optim2")
results_dir2 <- "data/pondTurt/pallida_sdm_optim2/"
#folder with pallida SDM
rast_dir <- "data/pondTurt/pall_sdm/"

## GA.prep and jl.prep


GA.inputs <- GA.prep(ASCII.dir = rast_dir, 
                     Results.dir = results_dir, parallel = 2)

GA.inputs2 <- GA.prep(ASCII.dir = rast_dir,
                      Results.dir = results_dir2, parallel = 2)

jl.inputs <- jl.prep(n.Pops = length(pall_pop), 
                     response = lower(pall_fst), 
                     CS_Point.File = pall_pop, JULIA_HOME = JULIA_HOME)


## optimize resistance surface

pallida_sdm_optim <- SS_optim(jl.inputs = jl.inputs, GA.inputs = GA.inputs)

pallida_sdm_optim2 <- SS_optim(jl.inputs = jl.inputs, GA.inputs = GA.inputs2)


## re-do original pallida ss optimize with 3 cores

dir.create("data/pondTurt/pallida_ssoptim_new")
results_dir <- "data/pondTurt/pallida_ssoptim_new/"

rast_dir <- "data/pondTurt/pall_surfaces/"


GA.inputs <- GA.prep(ASCII.dir = rast_dir, 
                     Results.dir = results_dir, parallel = 3)


jl.inputs <- jl.prep(n.Pops = length(pall_pop), 
                     response = lower(pall_chord), 
                     CS_Point.File = pall_pop, JULIA_HOME = JULIA_HOME)


## optimize resistance surface

pallida_ss_optim <- SS_optim(jl.inputs = jl.inputs, GA.inputs = GA.inputs)


## do MS_optim with all three

dir.create("data/pondTurt/pallida_msoptim_all")
results_dir <- "data/pondTurt/pallida_msoptim_all/"

rast_dir <- "data/pondTurt/pall_surfaces/"


GA.inputs <- GA.prep(ASCII.dir = rast_dir, 
                     Results.dir = results_dir, parallel = 2)


jl.inputs <- jl.prep(n.Pops = length(pall_pop), 
                     response = lower(pall_chord), 
                     CS_Point.File = pall_pop, JULIA_HOME = JULIA_HOME)




pallida_ms_optim <- MS_optim(jl.inputs = jl.inputs, GA.inputs = GA.inputs)
