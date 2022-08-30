#script for rasterize and other environmental data prep
# devtools::install_github("nspope/corMLPE")
# devtools::install_github("nspope/radish")

library(corMLPE)
library(radish)
library(sf)
library(tidyverse)
library(terra)

###bring in shapefiles and prep with crop and repro

#go to where the file is located.
#this can be done in read_sf call, but running into issues with current version
setwd("./env_data/NorWeST_PredictedStreamTempLines_WACoast_Aug/")
norwest <- read_sf(dsn = ".", layer = "NorWeST_PredictedStreamTempLines_WACoast_Aug")

setwd("../TerrainWorks_Dung_Elwha_Corrected/TerrainWorks_Dung_Elwha_Corrected/")
tworks <- read_sf(dsn = ".", layer = "Terrainworks_Dung_Elwha")

#from Amy is a template file for cropping both layers and set crs
setwd("../../../env_data/Elwha_streams/")
elwha <- read_sf(dsn = ".", layer = "Elwha_streams")

#note different crs
st_crs(elwha)
st_crs(norwest)
st_crs(tworks)

#plot of raw
ggplot() + 
  geom_sf(data = norwest, size = 1, color = "black") + 
  ggtitle("AOI Boundary Plot") + 
  coord_sf() +
  theme_bw()

ggplot() + 
  geom_sf(data = tworks, size = 1, color = "black") + 
  ggtitle("AOI Boundary Plot") + 
  coord_sf() +
  theme_bw()

#need to be reprojected to match amy's elwha file
tworks_repro <- st_transform(tworks, st_crs(elwha))
norwest_repro <- st_transform(norwest, st_crs(elwha))

#plot repro terrainworks with norwest to make sure repro worked
# ggplot() + 
#   geom_sf(data = tworks_repro, size = 1, color = "blue") + 
#   geom_sf(data = norwest, size = 1, color = "black") +
#   ggtitle("AOI Boundary Plot") + 
#   coord_sf() +
#   theme_bw()

#norwest needs to be cropped next
norwest_crop <- st_crop(norwest_repro, elwha)
tworks_crop <- st_crop(tworks_repro, elwha)

#check cropped norwest
ggplot() +
  geom_sf(data = tworks_crop, size = 1, color = "blue") +
  geom_sf(data = norwest_crop, size = 1, color = "black") +
  ggtitle("AOI Boundary Plot") +
  coord_sf() +
  theme_bw()


####next is create a raster from each of the values we are interested in

#make spatvector
norwest_crop <- vect(norwest_crop)

#make raster template
ra <- rast(norwest_crop)

#set resolution of raster template
res(ra) <- 1/1000

#and then rasterize
aug_st_temp <- terra::rasterize(norwest_crop, ra, field = "S2_02_11")

#df for ggplot
aug_st_temp_df <- as.data.frame(aug_st_temp, xy = TRUE)

#plot to check
ggplot() +
  geom_raster(data = aug_st_temp_df, aes(x = x, y = y, fill = S2_02_11)) +
  scale_fill_continuous(type = "viridis", na.value = "white") +
  ggtitle("S2_02_11") +
  coord_sf() +
  theme_bw()

#vector
tworks_crop <- vect(tworks_crop)

#and then rasterize
tworks_flow <- terra::rasterize(tworks_crop, ra, field = "FlowVel")

#df for ggplot
tworks_flow_df <- as.data.frame(tworks_flow, xy = TRUE)

#plot
ggplot() +
  geom_raster(data = tworks_flow_df, aes(x = x, y = y, fill = FlowVel)) +
  scale_fill_continuous(type = "viridis", na.value = "white") +
  ggtitle("S2_02_11") +
  coord_sf() +
  theme_bw()



#this all works, so now just need to pull in other pieces
#and then rasterize
tworks_flow <- terra::rasterize(tworks_crop, ra, field = "FlowVel")
tworks_flow_df <- as.data.frame(tworks_flow, xy = TRUE)
tworks_bfq <- terra::rasterize(tworks_crop, ra, field = "BFQ")
tworks_bfq_df <- as.data.frame(tworks_bfq, xy = TRUE)
aug_st_temp <- terra::rasterize(norwest_crop, ra, field = "S2_02_11")
aug_st_temp_df <- as.data.frame(aug_st_temp, xy = TRUE)
slope <- terra::rasterize(norwest_crop, ra, field = "SLOPE")
slope_df <- as.data.frame(slope, xy = TRUE)
canopy <- terra::rasterize(norwest_crop, ra, field = "CANOPY")
canopy_df <- as.data.frame(canopy, xy = TRUE)
cumdrainag <- terra::rasterize(norwest_crop, ra, field = "CUMDRAINAG")
cumdrainag_df <- as.data.frame(cumdrainag, xy = TRUE)

#convert to raster file type
tworks_flow <- raster::raster(tworks_flow)
tworks_bfq <- raster::raster(tworks_bfq)
aug_st_temp <- raster::raster(aug_st_temp)
slope <- raster::raster(slope)
canopy <- raster::raster(canopy)
cumdrainag <- raster::raster(cumdrainag)

#one other step is to remove -9999 and make NA for all NW temp data
aug_st_temp <- raster::reclassify(aug_st_temp, cbind(-Inf, 0, NA), right=FALSE)
slope <- raster::reclassify(slope, cbind(-Inf, 0, NA), right=FALSE)
canopy <- raster::reclassify(canopy, cbind(-Inf, 0, NA), right=FALSE)
cumdrainag <- raster::reclassify(cumdrainag, cbind(-Inf, 0, NA), right=FALSE)

#stack raster for radish
env_stack <- raster::stack(tworks_flow, tworks_bfq, aug_st_temp, slope, canopy, cumdrainag)

#dataframe from the stack for checking final alignment with plot
df1 <- as.data.frame(env_stack$FlowVel, xy = TRUE)
df2 <- as.data.frame(env_stack$S2_02_11, xy = TRUE)

#remove NA to speed up the plot
df1 <- df1 %>% drop_na()
df2 <- df2 %>% drop_na()

#plot the two.
#note that there are some weird missing data bits norwest around where the dam is
ggplot() + 
  geom_raster(data = df1, aes(x = x, y = y)) +
  geom_tile(data = df2, aes(x = x, y = y, colour = S2_02_11)) +
  coord_sf() +
  theme_bw()

#save raster stack as object
saveRDS(env_stack, file = "../../outputs/raster_env_stack.RDS")
