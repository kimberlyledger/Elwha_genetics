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
# Genetic data
cn_neutral_gd_pop <- readRDS("popdist_post_neutral_samplingsite_Chinook.rds")
cn_neutral_gd_pop <- as.matrix(cn_neutral_gd_pop)

# Load and aggregate raster data
covariates1 <- readRDS("raster_env_stack2.RDS")
covariates2 <- readRDS("riverscape_variables.RDS")
#covariates1_lowres <- terra::aggregate(covariates1, fact=5, na.rm = TRUE)
#covariates2_lowres <- terra::aggregate(covariates2, fact=5, na.rm = TRUE)
covariates3 <- readRDS("MWMT.RDS")
#covariates3_lowres <- terra::aggregate(covariates3, fact=5, na.rm = TRUE)

# Coordinates data
cn_sites_pop <- read.csv("Chinook_LongLat_pop.csv")
cn_sites_pop_ord <- cn_sites_pop %>% arrange(desc(Lat)) #arrange sites by descending latitude
cn_sites_coords_pop_raw <- SpatialPoints(coords = cn_sites_pop_ord[,2:3], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
#cn_sites_coords_pop <- spTransform(cn_sites_coords_pop_raw, crs(covariates1_lowres))
cn_sites_coords_pop <- spTransform(cn_sites_coords_pop_raw, crs(covariates1))

# Reorder gd matrices according to the descending latitude
cn_order <- cn_sites_pop_ord$Pop_ID
cn_neutral_gd_pop_ord <- cn_neutral_gd_pop[order(factor(rownames(cn_neutral_gd_pop), levels = cn_order)),order(factor(colnames(cn_neutral_gd_pop), levels = cn_order))] 

# Separate each raster layer for ResistanceGA functions later
# Assign "10 * maximum value of each layer" for NA values
#FlowVel_lr <- covariates1_lowres[[1]]
#FlowVel_lr@data@values[is.nan(FlowVel_lr@data@values)]<-NA
#FlowVel_lr[is.na(FlowVel_lr)] <- 10*maxValue(FlowVel_lr)
#writeRaster(FlowVel_lr, filename = paste0(cnibr.dir,"FlowVel_lr.asc"), overwrite = TRUE)
#writeRaster(FlowVel_lr, filename = paste0(shibr.dir,"FlowVel_lr.asc"), overwrite = TRUE)
#cov_flowvel <- stack(FlowVel_lr)

FlowVel <- covariates1[[1]]
FlowVel@data@values[is.nan(FlowVel@data@values)]<-NA
FlowVel[is.na(FlowVel)] <- 10*maxValue(FlowVel)
writeRaster(FlowVel, filename = paste0(cnibr.dir,"FlowVel.asc"), overwrite = TRUE)
writeRaster(FlowVel, filename = paste0(shibr.dir,"FlowVel.asc"), overwrite = TRUE)

#BFQ_lr <- covariates1_lowres[[2]]
#BFQ_lr@data@values[is.nan(BFQ_lr@data@values)]<-NA
#BFQ_lr[is.na(BFQ_lr)] <- 10*maxValue(BFQ_lr)
#writeRaster(BFQ_lr, filename = paste0(cnibr.dir,"BFQ_lr.asc"), overwrite = TRUE)
#writeRaster(BFQ_lr, filename = paste0(shibr.dir,"BFQ_lr.asc"), overwrite = TRUE)
#cov_bfq <- stack(BFQ_lr)

BFQ <- covariates1[[2]]
BFQ@data@values[is.nan(BFQ@data@values)]<-NA
BFQ[is.na(BFQ)] <- 10*maxValue(BFQ)
writeRaster(BFQ, filename = paste0(cnibr.dir,"BFQ.asc"), overwrite = TRUE)
writeRaster(BFQ, filename = paste0(shibr.dir,"BFQ.asc"), overwrite = TRUE)

#SLOPE_lr <- covariates1_lowres[[3]]
#SLOPE_lr@data@values[is.nan(SLOPE_lr@data@values)]<-NA
#SLOPE_lr[is.na(SLOPE_lr)] <- 10*maxValue(SLOPE_lr)
#writeRaster(SLOPE_lr, filename = paste0(cnibr.dir,"SLOPE_lr.asc"), overwrite = TRUE)
#writeRaster(SLOPE_lr, filename = paste0(shibr.dir,"SLOPE_lr.asc"), overwrite = TRUE)
#cov_slope <- stack(SLOPE_lr)

SLOPE <- covariates1[[3]]
SLOPE@data@values[is.nan(SLOPE@data@values)]<-NA
SLOPE[is.na(SLOPE)] <- 10*maxValue(SLOPE)
writeRaster(SLOPE, filename = paste0(cnibr.dir,"SLOPE.asc"), overwrite = TRUE)
writeRaster(SLOPE, filename = paste0(shibr.dir,"SLOPE.asc"), overwrite = TRUE)

#PRECIP_lr <- covariates1_lowres[[4]]
#PRECIP_lr@data@values[is.nan(PRECIP_lr@data@values)]<-NA
#PRECIP_lr[is.na(PRECIP_lr)] <- 10*maxValue(PRECIP_lr)
#writeRaster(PRECIP_lr, filename = paste0(cnibr.dir,"PRECIP_lr.asc"), overwrite = TRUE)
#writeRaster(PRECIP_lr, filename = paste0(shibr.dir,"PRECIP_lr.asc"), overwrite = TRUE)
#cov_precip <- stack(PRECIP_lr)

PRECIP <- covariates1[[4]]
PRECIP@data@values[is.nan(PRECIP@data@values)]<-NA
PRECIP[is.na(PRECIP)] <- 10*maxValue(PRECIP)
writeRaster(PRECIP, filename = paste0(cnibr.dir,"PRECIP.asc"), overwrite = TRUE)
writeRaster(PRECIP, filename = paste0(shibr.dir,"PRECIP.asc"), overwrite = TRUE)

#CANOPY_lr <- covariates1_lowres[[5]]
#CANOPY_lr@data@values[is.nan(CANOPY_lr@data@values)]<-NA
#CANOPY_lr[is.na(CANOPY_lr)] <- 10*maxValue(CANOPY_lr)
#writeRaster(CANOPY_lr, filename = paste0(cnibr.dir,"CANOPY_lr.asc"), overwrite = TRUE)
#writeRaster(CANOPY_lr, filename = paste0(shibr.dir,"CANOPY_lr.asc"), overwrite = TRUE)
#cov_canopy <- stack(CANOPY_lr)

CANOPY <- covariates1[[5]]
CANOPY@data@values[is.nan(CANOPY@data@values)]<-NA
CANOPY[is.na(CANOPY)] <- 10*maxValue(CANOPY)
writeRaster(CANOPY, filename = paste0(cnibr.dir,"CANOPY.asc"), overwrite = TRUE)
writeRaster(CANOPY, filename = paste0(shibr.dir,"CANOPY.asc"), overwrite = TRUE)

#IP_Chinook_lr <- covariates1_lowres[[6]]
#IP_Chinook_lr@data@values[is.nan(IP_Chinook_lr@data@values)]<-NA
#IP_Chinook_lr[is.na(IP_Chinook_lr)] <- 10*maxValue(IP_Chinook_lr)
#writeRaster(IP_Chinook_lr, filename = paste0(cnibr.dir,"IP_Chinook_lr.asc"), overwrite = TRUE)
#cov_ipchinook <- stack(IP_Chinook_lr)

IP_Chinook <- covariates1[[6]]
IP_Chinook@data@values[is.nan(IP_Chinook@data@values)]<-NA
IP_Chinook[is.na(IP_Chinook)] <- 10*maxValue(IP_Chinook)
writeRaster(IP_Chinook, filename = paste0(cnibr.dir,"IP_Chinook.asc"), overwrite = TRUE)

#IP_Steelhd_lr <- covariates1_lowres[[7]]
#IP_Steelhd_lr@data@values[is.nan(IP_Steelhd_lr@data@values)]<-NA
#IP_Steelhd_lr[is.na(IP_Steelhd_lr)] <- 10*maxValue(IP_Steelhd_lr)
#writeRaster(IP_Steelhd_lr, filename = paste0(shibr.dir,"IP_Steelhd_lr.asc"), overwrite = TRUE)
#cov_ipsteelhd <- stack(IP_Steelhd_lr)

IP_Steelhd <- covariates1[[7]]
IP_Steelhd@data@values[is.nan(IP_Steelhd@data@values)]<-NA
IP_Steelhd[is.na(IP_Steelhd)] <- 10*maxValue(IP_Steelhd)
writeRaster(IP_Steelhd, filename = paste0(shibr.dir,"IP_Steelhd.asc"), overwrite = TRUE)

# Define regions to fill in "0" values in Pool_freq_lr and logjams_10_lr based on spawnable_lr
#spawnable_lr <- covariates2_lowres[[4]]
#spawnable_lr@data@values[is.nan(spawnable_lr@data@values)]<-NA

spawnable <- covariates2[[4]]
spawnable@data@values[is.nan(spawnable@data@values)]<-NA

#Pool_freq_lr <- covariates2_lowres[[1]]
#Pool_freq_lr@data@values[is.nan(Pool_freq_lr@data@values)]<-NA
#Pool_freq_lr[is.na(Pool_freq_lr) & !is.na(spawnable_lr)] <- 0
#Pool_freq_lr[is.na(Pool_freq_lr)] <- 10*maxValue(Pool_freq_lr)
#writeRaster(Pool_freq_lr, filename = paste0(cnibr.dir,"Pool_freq_lr.asc"), overwrite = TRUE)
#writeRaster(Pool_freq_lr, filename = paste0(shibr.dir,"Pool_freq_lr.asc"), overwrite = TRUE)
#cov_poolfreq <- stack(Pool_freq_lr)

Pool_freq <- covariates2[[1]]
Pool_freq@data@values[is.nan(Pool_freq@data@values)]<-NA
Pool_freq[is.na(Pool_freq) & !is.na(spawnable)] <- 0
Pool_freq[is.na(Pool_freq)] <- 10*maxValue(Pool_freq)
writeRaster(Pool_freq, filename = paste0(cnibr.dir,"Pool_freq.asc"), overwrite = TRUE)
writeRaster(Pool_freq, filename = paste0(shibr.dir,"Pool_freq.asc"), overwrite = TRUE)

#logjams_10_lr <- covariates2_lowres[[2]]
#logjams_10_lr@data@values[is.nan(logjams_10_lr@data@values)]<-NA
#logjams_10_lr[is.na(logjams_10_lr) & !is.na(spawnable_lr)] <- 0
#logjams_10_lr[is.na(logjams_10_lr)] <- 10*maxValue(logjams_10_lr)
#writeRaster(logjams_10_lr, filename = paste0(cnibr.dir,"logjams_10_lr.asc"), overwrite = TRUE)
#writeRaster(logjams_10_lr, filename = paste0(shibr.dir,"logjams_10_lr.asc"), overwrite = TRUE)
#cov_logjams10 <- stack(logjams_10_lr)

logjams_10 <- covariates2[[2]]
logjams_10@data@values[is.nan(logjams_10@data@values)]<-NA
logjams_10[is.na(logjams_10) & !is.na(spawnable)] <- 0
logjams_10[is.na(logjams_10)] <- 10*maxValue(logjams_10)
writeRaster(logjams_10, filename = paste0(cnibr.dir,"logjams_10.asc"), overwrite = TRUE)
writeRaster(logjams_10, filename = paste0(shibr.dir,"logjams_10.asc"), overwrite = TRUE)

#spawnable_lr[is.na(spawnable_lr)] <- 10*maxValue(spawnable_lr)
#writeRaster(spawnable_lr, filename = paste0(cnibr.dir,"spawnable_lr.asc"), overwrite = TRUE)
#writeRaster(spawnable_lr, filename = paste0(shibr.dir,"spawnable_lr.asc"), overwrite = TRUE)
#cov_spawnable <- stack(spawnable_lr)

spawnable[is.na(spawnable)] <- 10*maxValue(spawnable)
writeRaster(spawnable, filename = paste0(cnibr.dir,"spawnable.asc"), overwrite = TRUE)
writeRaster(spawnable, filename = paste0(shibr.dir,"spawnable.asc"), overwrite = TRUE)

#MWAT_lr <- covariates3_lowres[[1]]
#MWAT_lr@data@values[is.nan(MWAT_lr@data@values)]<-NA
#MWAT_lr[is.na(MWAT_lr)] <- 10*maxValue(MWAT_lr)
#writeRaster(MWAT_lr, filename = paste0(cnibr.dir,"MWAT_lr.asc"), overwrite = TRUE)
#writeRaster(MWAT_lr, filename = paste0(shibr.dir,"MWAT_lr.asc"), overwrite = TRUE)
#cov_mwat <- stack(MWAT_lr)

MWMT <- covariates3[[1]]
MWMT@data@values[is.nan(MWMT@data@values)]<-NA
MWMT[is.na(MWMT)] <- 10*maxValue(MWMT)
writeRaster(MWMT, filename = paste0(cnibr.dir,"MWMT.asc"), overwrite = TRUE)
writeRaster(MWMT, filename = paste0(shibr.dir,"MWMT.asc"), overwrite = TRUE)

# Make raster stacks for downstream analyses
#covariates_lowres_ibr_cn <- stack(FlowVel_lr, IP_Chinook_lr, BFQ_lr, MWAT_lr, PRECIP_lr, SLOPE_lr, CANOPY_lr, Pool_freq_lr, logjams_10_lr, spawnable_lr) #full model covariates for Chinook
#covariates_lowres_ibrsb_cn <- stack(FlowVel_lr, IP_Chinook_lr, BFQ_lr, MWAT_lr, PRECIP_lr, SLOPE_lr, CANOPY_lr) #subset model covariates excluding "riverscape variables" which has limited range

covariates_ibr_cn <- stack(FlowVel, IP_Chinook, MWMT, PRECIP, SLOPE, CANOPY, Pool_freq, logjams_10, spawnable) #full model covariates for Chinook
covariates_ibrsb_cn <- stack(FlowVel, IP_Chinook, MWMT, PRECIP, SLOPE, CANOPY) #subset model covariates excluding "riverscape variables" which has limited range

#covariates_lowres_flow <- stack(FlowVel_lr, BFQ_lr) #flow model covariates
#covariates_lowres_clim <- stack(MWAT_lr, PRECIP_lr) #climate model covariates
#covariates_lowres_cnhbt <- stack(IP_Chinook_lr, SLOPE_lr, CANOPY_lr, Pool_freq_lr, logjams_10_lr, spawnable_lr) #chinook habitat model covariates
#covariates_lowres_cnhbtsb <- stack(IP_Chinook_lr, SLOPE_lr, CANOPY_lr) #chinook habitat subset model covariates 

covariates_flow <- stack(FlowVel) #flow model covariates
covariates_clim <- stack(MWMT, PRECIP) #climate model covariates
covariates_cnhbt <- stack(IP_Chinook, SLOPE, CANOPY, Pool_freq, logjams_10, spawnable) #chinook habitat model covariates
covariates_cnhbtsb <- stack(IP_Chinook, SLOPE, CANOPY) #chinook habitat subset model covariates 

#OUT_DIST_lr <- covariates1_lowres[[8]]
#OUT_DIST_lr@data@values[is.nan(OUT_DIST_lr@data@values)]<-NA
#OUT_DIST_lr[is.na(OUT_DIST_lr)] <- 10*maxValue(OUT_DIST_lr)
#writeRaster(OUT_DIST_lr, filename = paste0(ibd.dir,"OUT_DIST_lr.asc"), overwrite = TRUE)

OUT_DIST <- covariates1[[8]]
OUT_DIST@data@values[is.nan(OUT_DIST@data@values)]<-NA
OUT_DIST[is.na(OUT_DIST)] <- 10*maxValue(OUT_DIST)
writeRaster(OUT_DIST, filename = paste0(ibd.dir,"OUT_DIST.asc"), overwrite = TRUE)

#covariates_lowres_ibd <- stack(OUT_DIST_lr)

covariates_ibd <- stack(OUT_DIST)


# Move point occurrences falling in raster cells without data (i.e. NA) to the nearest raster cell with data
# They might not move because all rater cells have values (at least "zero") from the above.
#cn_moved_coords_pop <- points2nearestcell(
#  locs = cn_sites_coords_pop,
#  ras = covariates_lowres_ibr_cn,
#  layer = 10,
#  move = TRUE,
#  distance = NULL,
#  showchanges = TRUE,
#  showmap = TRUE,
#  leaflet = FALSE
#)
#write.table(cn_moved_coords_pop,file=paste0(write.dir,"cn_moved_coords_pop.txt"),sep="\t",col.names=F,row.names=F)
#cn.locales <- read.table("/scratch/bell/jeon96/LGC/IBR/rga/cn_moved_coords_pop.txt")
#cn.locales <- SpatialPoints(cn.locales[,c(1,2)])
#crs(cn.locales) <- crs(covariates_lowres_ibr_cn$FlowVel) #define projection

cn_moved_coords_pop <- points2nearestcell(
  locs = cn_sites_coords_pop,
  ras = covariates_ibr_cn,
  layer = 10,
  move = TRUE,
  distance = NULL,
  showchanges = TRUE,
  showmap = TRUE,
  leaflet = FALSE
)
write.table(cn_moved_coords_pop,file=paste0(write.dir,"cn_moved_coords_pop.txt"),sep="\t",col.names=F,row.names=F)
cn.locales <- read.table("/scratch/bell/jeon96/LGC/IBR/rga/cn_moved_coords_pop.txt")
cn.locales <- SpatialPoints(cn.locales[,c(1,2)])
crs(cn.locales) <- crs(covariates_ibr_cn$FlowVel) #define projection

print("Data loading and preprocessing has been finished.")

## Chinook_neutral
# Change paths below to yours
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/ibr/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/ibrsb/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/flow/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/clim/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/hbt/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/hbtsb/"))
dir.create(file.path("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/ibd/"))
cnnt_ibr.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/cnnt/ibr/"
cnnt_ibrsb.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/cnnt/ibrsb/"
cnnt_flow.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/cnnt/flow/"
cnnt_clim.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/cnnt/clim/"
cnnt_hbt.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/cnnt/hbt/"
cnnt_hbtsb.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/cnnt/hbtsb/"
cnnt_ibd.dir <- "/scratch/bell/jeon96/LGC/IBR/rga/cnnt/ibd/"

cnnt_gdist.inputs <- gdist.prep(n.Pops = length(cn.locales),
                                samples = cn.locales,
                                response = as.vector(lower(cn_neutral_gd_pop_ord)),
                                method = 'commuteDistance')
                                
# Chinook neutral IBR with 1 replicate, 1 iteration (i.e., no bootstrapping)
# individual and full model
#cnnt_ibr_GA.inputs <- GA.prep(method = "LL",
#                              ASCII.dir = covariates_lowres_ibr_cn,
#                              scale = FALSE,
#                              Results.dir = cnnt_ibr.dir,
#                              parallel = 64)

cnnt_ibr_GA.inputs <- GA.prep(method = "LL",
                              ASCII.dir = covariates_ibr_cn,
                              scale = FALSE,
                              Results.dir = cnnt_ibr.dir,
                              parallel = 64)
                              
#cnnt_ibr_SS_results <- SS_optim(gdist.inputs = cnnt_gdist.inputs,
#                                GA.inputs = cnnt_ibr_GA.inputs,
#                                diagnostic_plots = FALSE)
#saveRDS(cnnt_ibr_SS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_ibr_SS_results.RDS")                                
                                                           
cnnt_ibr_MS_results <- MS_optim(gdist.inputs = cnnt_gdist.inputs,
                                GA.inputs = cnnt_ibr_GA.inputs,
                                diagnostic_plots = FALSE)
saveRDS(cnnt_ibr_MS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_ibr_MS_results.RDS")     
print("Full ibr model has been optimized.")

# full model (subset - excluding layers from "riverscape variables")
#cnnt_ibrsb_GA.inputs <- GA.prep(method = "LL",
#                                ASCII.dir = covariates_lowres_ibrsb_cn,
#                                scale = FALSE,
#                                Results.dir = cnnt_ibrsb.dir,
#                                parallel = 64)

cnnt_ibrsb_GA.inputs <- GA.prep(method = "LL",
                                ASCII.dir = covariates_ibrsb_cn,
                                scale = FALSE,
                                Results.dir = cnnt_ibrsb.dir,
                                parallel = 64)
                                                                
cnnt_ibrsb_MS_results <- MS_optim(gdist.inputs = cnnt_gdist.inputs,
                                  GA.inputs = cnnt_ibrsb_GA.inputs,
                                  diagnostic_plots = FALSE)
saveRDS(cnnt_ibrsb_MS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_ibrsb_MS_results.RDS")     
print("Subset ibr model has been optimized.")
                                
# flow model                            
#cnnt_flow_GA.inputs <- GA.prep(method = "LL",
#                               ASCII.dir = covariates_lowres_flow,
#                               scale = FALSE,
#                               Results.dir = cnnt_flow.dir,
#                               parallel = 64)

cnnt_flow_GA.inputs <- GA.prep(method = "LL",
                               ASCII.dir = covariates_flow,
                               scale = FALSE,
                               Results.dir = cnnt_flow.dir,
                               parallel = 64)
                               
cnnt_flow_SS_results <- SS_optim(gdist.inputs = cnnt_gdist.inputs,
                                 GA.inputs = cnnt_flow_GA.inputs,
                                 diagnostic_plots = FALSE)
saveRDS(cnnt_flow_SS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_flow_SS_results.RDS")    
print("Flow model has been optimized.")
                             
# clim model                            
#cnnt_clim_GA.inputs <- GA.prep(method = "LL",
#                               ASCII.dir = covariates_lowres_clim,
#                               scale = FALSE,
#                               Results.dir = cnnt_clim.dir,
#                               parallel = 64)

cnnt_clim_GA.inputs <- GA.prep(method = "LL",
                               ASCII.dir = covariates_clim,
                               scale = FALSE,
                               Results.dir = cnnt_clim.dir,
                               parallel = 64)
                               
cnnt_clim_MS_results <- MS_optim(gdist.inputs = cnnt_gdist.inputs,
                                 GA.inputs = cnnt_clim_GA.inputs,
                                 diagnostic_plots = FALSE)
saveRDS(cnnt_clim_MS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_clim_MS_results.RDS")
print("Climate model has been optimized.")
                                 
# habitat model (full)
#cnnt_hbt_GA.inputs <- GA.prep(method = "LL",
#                              ASCII.dir = covariates_lowres_cnhbt,
#                              scale = FALSE,
#                              Results.dir = cnnt_hbt.dir,
#                              parallel = 64)

cnnt_hbt_GA.inputs <- GA.prep(method = "LL",
                              ASCII.dir = covariates_cnhbt,
                              scale = FALSE,
                              Results.dir = cnnt_hbt.dir,
                              parallel = 64)
                              
cnnt_hbt_MS_results <- MS_optim(gdist.inputs = cnnt_gdist.inputs,
                                GA.inputs = cnnt_hbt_GA.inputs,
                                diagnostic_plots = FALSE)
saveRDS(cnnt_hbt_MS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_hbt_MS_results.RDS")
print("Habitat model has been optimized.")
                                 
# habitat model (subset)
#cnnt_hbtsb_GA.inputs <- GA.prep(method = "LL",
#                                ASCII.dir = covariates_lowres_cnhbtsb,
#                                scale = FALSE,
#                                Results.dir = cnnt_hbtsb.dir,
#                                parallel = 64)

cnnt_hbtsb_GA.inputs <- GA.prep(method = "LL",
                                ASCII.dir = covariates_cnhbtsb,
                                scale = FALSE,
                                Results.dir = cnnt_hbtsb.dir,
                                parallel = 64)
                                
cnnt_hbtsb_MS_results <- MS_optim(gdist.inputs = cnnt_gdist.inputs,
                                  GA.inputs = cnnt_hbtsb_GA.inputs,
                                  diagnostic_plots = FALSE)
saveRDS(cnnt_hbtsb_MS_results, file="/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_hbtsb_MS_results.RDS")
print("Subset habitat model has been optimized.")
                                 
# Chinook neutral IBD with 1 replicate, 1 iteration (i.e., no bootstrapping)
#cnnt_ibd_GA.inputs <- GA.prep(method = "LL",
#                              ASCII.dir = covariates_lowres_ibd,
#                              scale = FALSE,
#                              Results.dir = cnnt_ibd.dir,
#                              parallel = 64)                              

cnnt_ibd_GA.inputs <- GA.prep(method = "LL",
                              ASCII.dir = covariates_ibd,
                              scale = FALSE,
                              Results.dir = cnnt_ibd.dir,
                              parallel = 64) 
                              
cnnt_ibd_results <- SS_optim(gdist.inputs = cnnt_gdist.inputs,
                             GA.inputs = cnnt_ibd_GA.inputs,
                             diagnostic_plots = FALSE) 
saveRDS(cnnt_ibd_results, file="/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_ibd_results.RDS")
print("IBD model has been optimized.")


# Boostrapping
# Extract relevant components from optimization outputs
# Make a list of cost/resistance distance matrices
cnnt_ibr_MS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_ibr_MS_results.RDS")
cnnt_ibrsb_MS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_ibrsb_MS_results.RDS")
cnnt_flow_SS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_flow_SS_results.RDS")
cnnt_clim_MS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_clim_MS_results.RDS")
cnnt_hbt_MS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_hbt_MS_results.RDS")
cnnt_hbtsb_MS_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_hbtsb_MS_results.RDS")
cnnt_ibd_results <- readRDS("/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_ibd_results.RDS")

cnnt_mat.list <- c(cnnt_ibr_MS_results$cd,
                   cnnt_ibrsb_MS_results$cd,
                   cnnt_flow_SS_results$cd,
                   cnnt_clim_MS_results$cd,
                   cnnt_hbt_MS_results$cd,
                   cnnt_hbtsb_MS_results$cd,
                   cnnt_ibd_results$cd)

cnnt_k <- rbind(cnnt_ibr_MS_results$k,
                cnnt_ibrsb_MS_results$k,
                cnnt_flow_SS_results$k,
                cnnt_clim_MS_results$k,
                cnnt_hbt_MS_results$k,
                cnnt_hbtsb_MS_results$k,
                cnnt_ibd_results$k)

# Create square distance matrix for response for use with
# the bootstrap function        
cnnt_response <- matrix(0, 10, 10) #number of populations
cnnt_response[lower.tri(cnnt_response)] <- lower(cn_neutral_gd_pop_ord)

# Run bootstrap
(cnnt_AIC.boot <- Resist.boot(mod.names = names(cnnt_mat.list),
                              dist.mat = cnnt_mat.list,
                              n.parameters = cnnt_k[,2],
                              sample.prop = 0.75,
                              iters = 1000,
                              obs = 10,
                              genetic.mat = cnnt_response))
saveRDS(cnnt_AIC.boot, file="/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_AIC_result.RDS")
write.table(cnnt_AIC.boot, file = "/scratch/bell/jeon96/LGC/IBR/rga/cnnt/cnnt_AIC_result.txt")
print("Bootstrapping has been finished.")
print("All RGA anayses has been finished for Chinook neutral genetic distance response.")