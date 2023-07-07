##Load ResistanceGA
library(ResistanceGA, lib.loc = "/home/jeon96/R/bell/4.1.2-gcc-9.3.0-rw7vp7m/")
library(raster, lib.loc = "/home/jeon96/R/bell/4.1.2-gcc-9.3.0-rw7vp7m/")
library(sf, lib.loc = "/home/jeon96/R/bell/4.1.2-gcc-9.3.0-rw7vp7m/")
library(rSDM, lib.loc = "/home/jeon96/R/bell/4.1.2-gcc-9.3.0-rw7vp7m/")
library(terra, lib.loc = "/home/jeon96/R/bell/4.1.2-gcc-9.3.0-rw7vp7m/")
library(dplyr)
library(parallel)
detectCores()

setwd(dir = "/scratch/bell/jeon96/LGC/IBR/")
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/ibr/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/ibd/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/cnibr/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/shibr/"))
write.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/"
cnibr.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/cnibr/"
shibr.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/shibr/"
ibd.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/ibd/"

##Load data
# Genetic data (steelhead adults)
shjuv_neutral_gd_pop <- readRDS("gen_steel_post_neutral_juvenile.rds")
shjuv_neutral_gd_pop <- as.matrix(shjuv_neutral_gd_pop)

# Load and aggregate raster data
covariates1 <- readRDS("raster_env_stack2.RDS")
covariates2 <- readRDS("riverscape_variables.RDS")
#covariates1_lowres <- terra::aggregate(covariates1, fact=5, na.rm = TRUE)
#covariates2_lowres <- terra::aggregate(covariates2, fact=5, na.rm = TRUE)
covariates3 <- readRDS("MWMT.RDS")
#covariates3_lowres <- terra::aggregate(covariates3, fact=5, na.rm = TRUE)

# Coordinates data
shjuv_sites_pop <- read.csv("Steelhead_juv_LongLat_pop2.csv")
shjuv_sites_pop_ord <- shjuv_sites_pop %>% arrange(desc(Lat)) #arrange sites by descending latitude
shjuv_sites_coords_pop_raw <- SpatialPoints(coords = shjuv_sites_pop_ord[,2:3], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
#sh_sites_coords_pop <- spTransform(sh_sites_coords_pop_raw, crs(covariates1_lowres))
shjuv_sites_coords_pop <- spTransform(shjuv_sites_coords_pop_raw, crs(covariates1))

# Reorder gd matrices according to the descending latitude
shjuv_order <- shjuv_sites_pop_ord$Pop_ID
shjuv_neutral_gd_pop_ord <- shjuv_neutral_gd_pop[order(factor(rownames(shjuv_neutral_gd_pop), levels = shjuv_order)),order(factor(colnames(shjuv_neutral_gd_pop), levels = shjuv_order))] 

# Separate each raster layer for ResistanceGA functions later
# Assign "10 * maximum value of each layer" for NA values
FlowVel <- covariates1[[1]]
FlowVel@data@values[is.nan(FlowVel@data@values)]<-NA
FlowVel[is.na(FlowVel)] <- 10*maxValue(FlowVel)
writeRaster(FlowVel, filename = paste0(cnibr.dir,"FlowVel.asc"), overwrite = TRUE)
writeRaster(FlowVel, filename = paste0(shibr.dir,"FlowVel.asc"), overwrite = TRUE)

BFQ <- covariates1[[2]]
BFQ@data@values[is.nan(BFQ@data@values)]<-NA
BFQ[is.na(BFQ)] <- 10*maxValue(BFQ)
writeRaster(BFQ, filename = paste0(cnibr.dir,"BFQ.asc"), overwrite = TRUE)
writeRaster(BFQ, filename = paste0(shibr.dir,"BFQ.asc"), overwrite = TRUE)

SLOPE <- covariates1[[3]]
SLOPE@data@values[is.nan(SLOPE@data@values)]<-NA
SLOPE[is.na(SLOPE)] <- 10*maxValue(SLOPE)
writeRaster(SLOPE, filename = paste0(cnibr.dir,"SLOPE.asc"), overwrite = TRUE)
writeRaster(SLOPE, filename = paste0(shibr.dir,"SLOPE.asc"), overwrite = TRUE)

PRECIP <- covariates1[[4]]
PRECIP@data@values[is.nan(PRECIP@data@values)]<-NA
PRECIP[is.na(PRECIP)] <- 10*maxValue(PRECIP)
writeRaster(PRECIP, filename = paste0(cnibr.dir,"PRECIP.asc"), overwrite = TRUE)
writeRaster(PRECIP, filename = paste0(shibr.dir,"PRECIP.asc"), overwrite = TRUE)

CANOPY <- covariates1[[5]]
CANOPY@data@values[is.nan(CANOPY@data@values)]<-NA
CANOPY[is.na(CANOPY)] <- 10*maxValue(CANOPY)
writeRaster(CANOPY, filename = paste0(cnibr.dir,"CANOPY.asc"), overwrite = TRUE)
writeRaster(CANOPY, filename = paste0(shibr.dir,"CANOPY.asc"), overwrite = TRUE)

IP_Chinook <- covariates1[[6]]
IP_Chinook@data@values[is.nan(IP_Chinook@data@values)]<-NA
IP_Chinook[is.na(IP_Chinook)] <- 10*maxValue(IP_Chinook)
writeRaster(IP_Chinook, filename = paste0(cnibr.dir,"IP_Chinook.asc"), overwrite = TRUE)

IP_Steelhd <- covariates1[[7]]
IP_Steelhd@data@values[is.nan(IP_Steelhd@data@values)]<-NA
IP_Steelhd[is.na(IP_Steelhd)] <- 10*maxValue(IP_Steelhd)
writeRaster(IP_Steelhd, filename = paste0(shibr.dir,"IP_Steelhd.asc"), overwrite = TRUE)

# Define regions to fill in "0" values in Pool_freq_lr and logjams_10_lr based on spawnable_lr
spawnable <- covariates2[[4]]
spawnable@data@values[is.nan(spawnable@data@values)]<-NA

Pool_freq <- covariates2[[1]]
Pool_freq@data@values[is.nan(Pool_freq@data@values)]<-NA
Pool_freq[is.na(Pool_freq) & !is.na(spawnable)] <- 0
Pool_freq[is.na(Pool_freq)] <- 10*maxValue(Pool_freq)
writeRaster(Pool_freq, filename = paste0(cnibr.dir,"Pool_freq.asc"), overwrite = TRUE)
writeRaster(Pool_freq, filename = paste0(shibr.dir,"Pool_freq.asc"), overwrite = TRUE)

logjams_10 <- covariates2[[2]]
logjams_10@data@values[is.nan(logjams_10@data@values)]<-NA
logjams_10[is.na(logjams_10) & !is.na(spawnable)] <- 0
logjams_10[is.na(logjams_10)] <- 10*maxValue(logjams_10)
writeRaster(logjams_10, filename = paste0(cnibr.dir,"logjams_10.asc"), overwrite = TRUE)
writeRaster(logjams_10, filename = paste0(shibr.dir,"logjams_10.asc"), overwrite = TRUE)

spawnable[is.na(spawnable)] <- 10*maxValue(spawnable)
writeRaster(spawnable, filename = paste0(cnibr.dir,"spawnable.asc"), overwrite = TRUE)
writeRaster(spawnable, filename = paste0(shibr.dir,"spawnable.asc"), overwrite = TRUE)

MWMT <- covariates3[[1]]
MWMT@data@values[is.nan(MWMT@data@values)]<-NA
MWMT[is.na(MWMT)] <- 10*maxValue(MWMT)
writeRaster(MWAT, filename = paste0(cnibr.dir,"MWAT.asc"), overwrite = TRUE)
writeRaster(MWAT, filename = paste0(shibr.dir,"MWAT.asc"), overwrite = TRUE)

# Make raster stacks for downstream analyses
covariates_ibr_sh <- stack(FlowVel, IP_Steelhd, MWMT, PRECIP, SLOPE, CANOPY, Pool_freq, logjams_10, spawnable) #full model covariates for Steelhead
covariates_ibrsb_sh <- stack(FlowVel, IP_Steelhd, MWMT, PRECIP, SLOPE, CANOPY) #subset model covariates excluding "riverscape variables" which has limited range

covariates_flow <- stack(FlowVel) #flow model covariates
covariates_clim <- stack(MWMT, PRECIP) #climate model covariates
covariates_shhbt <- stack(IP_Steelhd, SLOPE, CANOPY, Pool_freq, logjams_10, spawnable) #steelhead habitat model covariates
covariates_shhbtsb <- stack(IP_Steelhd, SLOPE, CANOPY) #steelhead habitat subset model covariates

OUT_DIST <- covariates1[[8]]
OUT_DIST@data@values[is.nan(OUT_DIST@data@values)]<-NA
OUT_DIST[is.na(OUT_DIST)] <- 10*maxValue(OUT_DIST)
writeRaster(OUT_DIST, filename = paste0(ibd.dir,"OUT_DIST.asc"), overwrite = TRUE)

#covariates_lowres_ibd <- stack(OUT_DIST_lr)

covariates_ibd <- stack(OUT_DIST)

# Move point occurrences falling in raster cells without data (i.e. NA) to the nearest raster cell with data
# They might not move because all rater cells have values (at least "zero") from the above.

shjuv_moved_coords_pop <- points2nearestcell(
  locs = shjuv_sites_coords_pop,
  ras = covariates_ibr_sh,
  layer = 10,
  move = TRUE,
  distance = NULL,
  showchanges = TRUE,
  showmap = TRUE,
  leaflet = FALSE
)
write.table(shjuv_moved_coords_pop,file=paste0(write.dir,"shjuv_moved_coords_pop.txt"),sep="\t",col.names=F,row.names=F)
shjuv.locales <- read.table("/scratch/bell/jeon96/LGC/IBR/rga/shjuv_moved_coords_pop.txt")
shjuv.locales <- SpatialPoints(shjuv.locales[,c(1,2)])
crs(shjuv.locales) <- crs(covariates_ibr_sh$FlowVel) #define projection

print("Data loading and preprocessing has been finished.")

## Steelhead_neutral_juvenile
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/ibr"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/ibrsb"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/flow"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/clim"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/hbt"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/hbtsb"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/ibd"))
shnt_juv_ibr.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/ibr/"
shnt_juv_ibrsb.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/ibrsb/"
shnt_juv_flow.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/flow/"
shnt_juv_clim.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/clim/"
shnt_juv_hbt.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/hbt/"
shnt_juv_hbtsb.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/hbtsb/"
shnt_juv_ibd.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/ibd/"

# Change the number of cores below
shnt_juv_gdist.inputs <- gdist.prep(n.Pops = length(shjuv.locales),
                                    samples = shjuv.locales,
                                    response = as.vector(lower(shjuv_neutral_gd_pop_ord)),
                                    method = 'commuteDistance')

# Change the number of cores in "parallel" below
# Steelhead juveniles neutral IBR with 1 replicate, 1 iteration (i.e., no bootstrapping)
# individual and full model
print("Full ibr model optimization started.")
shnt_juv_ibr_GA.inputs <- GA.prep(method = "LL",
                                  ASCII.dir = covariates_ibr_sh,
                                  scale = FALSE,
                                  Results.dir = shnt_juv_ibr.dir,
                                  parallel = 64)

#shnt_juv_ibr_SS_results <- SS_optim(gdist.inputs = shnt_juv_gdist.inputs,
#                                    GA.inputs = shnt_juv_ibr_GA.inputs,
#                                    diagnostic_plots = FALSE)
#saveRDS(shnt_juv_ibr_SS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_ibr_SS_results.RDS")
                                                           
shnt_juv_ibr_MS_results <- MS_optim(gdist.inputs = shnt_juv_gdist.inputs,
                                    GA.inputs = shnt_juv_ibr_GA.inputs,
                                    diagnostic_plots = FALSE)
saveRDS(shnt_juv_ibr_MS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_ibr_MS_results.RDS")
print("Full ibr model has been optimized.")

# full model (subset - excluding layers from "riverscape variables")
print("Subset ibr model optimization started.")
shnt_juv_ibrsb_GA.inputs <- GA.prep(method = "LL",
                                    ASCII.dir = covariates_ibrsb_sh,
                                    scale = FALSE,
                                    Results.dir = shnt_juv_ibrsb.dir,
                                    parallel = 64)
                                
shnt_juv_ibrsb_MS_results <- MS_optim(gdist.inputs = shnt_juv_gdist.inputs,
                                      GA.inputs = shnt_juv_ibrsb_GA.inputs,
                                      diagnostic_plots = FALSE)
saveRDS(shnt_juv_ibrsb_MS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_ibrsb_MS_results.RDS")
print("Subset ibr model has been optimized.")
                                
# flow model     
print("Flow model optimization started.")                        
shnt_juv_flow_GA.inputs <- GA.prep(method = "LL",
                                   ASCII.dir = covariates_flow,
                                   scale = FALSE,
                                   Results.dir = shnt_juv_flow.dir,
                                   parallel = 64)

shnt_juv_flow_SS_results <- SS_optim(gdist.inputs = shnt_juv_gdist.inputs,
                                     GA.inputs = shnt_juv_flow_GA.inputs,
                                     diagnostic_plots = FALSE)
saveRDS(shnt_juv_flow_SS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_flow_SS_results.RDS")
print("Flow model has been optimized.")
                                 
# clim model 
print("Climate model optimization started.")                        
shnt_juv_clim_GA.inputs <- GA.prep(method = "LL",
                                   ASCII.dir = covariates_clim,
                                   scale = FALSE,
                                   Results.dir = shnt_juv_clim.dir,
                                   parallel = 64)

shnt_juv_clim_MS_results <- MS_optim(gdist.inputs = shnt_juv_gdist.inputs,
                                     GA.inputs = shnt_juv_clim_GA.inputs,
                                     diagnostic_plots = FALSE)
saveRDS(shnt_juv_clim_MS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_clim_MS_results.RDS")
print("Climate model has been optimized.")
                                 
# habitat model (full)
print("Habitat model optimization started.")
shnt_juv_hbt_GA.inputs <- GA.prep(method = "LL",
                                  ASCII.dir = covariates_shhbt,
                                  scale = FALSE,
                                  Results.dir = shnt_juv_hbt.dir,
                                  parallel = 64)

shnt_juv_hbt_MS_results <- MS_optim(gdist.inputs = shnt_juv_gdist.inputs,
                                    GA.inputs = shnt_juv_hbt_GA.inputs,
                                    diagnostic_plots = FALSE)
saveRDS(shnt_juv_hbt_MS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_hbt_MS_results.RDS")
print("Habitat model has been optimized.")
                                 
# habitat model (subset)
print("Subset habitat model optimization started.")
shnt_juv_hbtsb_GA.inputs <- GA.prep(method = "LL",
                                    ASCII.dir = covariates_shhbtsb,
                                    scale = FALSE,
                                    Results.dir = shnt_juv_hbtsb.dir,
                                    parallel = 64)

shnt_juv_hbtsb_MS_results <- MS_optim(gdist.inputs = shnt_juv_gdist.inputs,
                                      GA.inputs = shnt_juv_hbtsb_GA.inputs,
                                      diagnostic_plots = FALSE)
saveRDS(shnt_juv_hbtsb_MS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_hbtsb_MS_results.RDS")
print("Subset habitat model has been optimized.")
                                 
# Steelhead juveniles neutral IBD with 1 replicate, 1 iteration (i.e., no bootstrapping)
print("IBD model optimization started.")
shnt_juv_ibd_GA.inputs <- GA.prep(method = "LL",
                                  ASCII.dir = covariates_ibd,
                                  scale = FALSE,
                                  Results.dir = shnt_juv_ibd.dir,
                                  parallel = 64)                              

shnt_juv_ibd_results <- SS_optim(gdist.inputs = shnt_juv_gdist.inputs,
                                 GA.inputs = shnt_juv_ibd_GA.inputs,
                                 diagnostic_plots = FALSE) 
saveRDS(shnt_juv_ibd_results, file="/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_ibd_results.RDS")
print("IBD model has been optimized.")


# Boostrapping
# Extract relevant components from optimization outputs
# Make a list of cost/resistance distance matrices
shnt_juv_ibr_MS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_ibr_MS_results.RDS")
shnt_juv_ibrsb_MS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_ibrsb_MS_results.RDS")
shnt_juv_flow_SS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_flow_SS_results.RDS")
shnt_juv_clim_MS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_clim_MS_results.RDS")
shnt_juv_hbt_MS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_hbt_MS_results.RDS")
shnt_juv_hbtsb_MS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_hbtsb_MS_results.RDS")
shnt_juv_ibd_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_ibd_results.RDS")

shnt_juv_mat.list <- c(shnt_juv_ibr_MS_results$cd,
                       shnt_juv_ibrsb_MS_results$cd,
                       shnt_juv_flow_SS_results$cd,
                       shnt_juv_clim_MS_results$cd,
                       shnt_juv_hbt_MS_results$cd,
                       shnt_juv_hbtsb_MS_results$cd,
                       shnt_juv_ibd_results$cd)

shnt_juv_k <- rbind(shnt_juv_ibr_MS_results$k,
                    shnt_juv_ibrsb_MS_results$k,
                    shnt_juv_flow_SS_results$k,
                    shnt_juv_clim_MS_results$k,
                    shnt_juv_hbt_MS_results$k,
                    shnt_juv_hbtsb_MS_results$k,
                    shnt_juv_ibd_results$k)

# Create square distance matrix for response for use with
# the bootstrap function        
shnt_juv_response <- matrix(0, 4, 4) #number of populations
shnt_juv_response[lower.tri(shnt_juv_response)] <- lower(shjuv_neutral_gd_pop_ord)

# Run bootstrap
(shnt_juv_AIC.boot <- Resist.boot(mod.names = names(shnt_juv_mat.list),
                                  dist.mat = shnt_juv_mat.list,
                                  n.parameters = shnt_juv_k[,2],
                                  sample.prop = 0.75,
                                  iters = 1000,
                                  obs = 4,
                                  genetic.mat = shnt_juv_response))
saveRDS(shnt_juv_AIC.boot, file = "/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_AIC_result.RDS")
write.table(shnt_juv_AIC.boot, file = "/scratch/bell/jeon96/LGC/IBR/rga/shnt_juv/shnt_juv_AIC_result.txt")
print("Bootstrapping has been finished.")
print("All RGA anayses has been finished for Steelhead/rainbow trout juvenile neutral genetic distance response.")