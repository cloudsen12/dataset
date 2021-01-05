#' Create points in cloudsen12
#' @author csaybar
#'
#' Script used to create manually points in cloudsen12
#'


library(tidyverse)
library(mapedit)
library(mapview)
library(raster)
library(tmap)
library(rgee)
library(sf)

source("src/utils.R")
data("World")
ee_Initialize("csaybar", gcs = TRUE)

# 1. metadata viz GLOBAL
land_use_metadata <- tibble(
  class = c("Barren", "Tropical Forest", "Temperated Forest",
            "Grass/Crop", "shrubland", "Snow/Ice/Lichen",
            "Urban", "Water", "Wetlands", "Ocean"),
  color = c("B4B4B4", "8DB400", "A0DC00", "F096FF",
            "FFBB22", "00E1FF", "FA0000", "0032C8",
            "0096A0", "0000FF"),
  value = 1:10
)
landuse_world <- ee$Image("users/csaybar/cloudsen2/world_landuse")
eo_compass <- ee$FeatureCollection("users/csaybar/cloudsen2/eo_compass")
eo_compass_img <- ee$Image()$float()$paint(eo_compass, "nmbrfsc")
palette <- viridisLite::viridis(10)

# previous points -------------------------------
pot_points_sf <- read_sf("data/cloudsen2_potential_points.geojson")

map_base <- Map$addLayer(landuse_world, list(min = 0, max = 10, palette = land_use_metadata$color), shown = FALSE) +
  Map$addLayer(eo_compass_img, list(min = 60, max = 100, palette = palette), legend = TRUE, name = "eocompass", opacity = 0.5)

map_base$rgee$tokens
# 2. Select Potential points by class ------------------------------------------
barren <- landuse_world$eq(1)

barren_map <- Map$addLayer(barren$updateMask(barren), list(min = 1, max = 10), name = "baren") +
  map_base
barren_map <- mapview(pot_points_sf, barren_map, zcol = "good", cex = 3, legend = FALSE)
barren_points <- mapedit::editMap(barren_map)


cloudsen2 <- barren_points$drawn['geometry']
cloudsen2$type <- "baren"



## Tropical Forest
tforest <- landuse_world$eq(2)
map_tforest <- Map$addLayer(tforest$updateMask(tforest), list(min = 1, max = 10), name = "tf") +
  map_base
map_tforest <- mapview(cloudsen2, map_tforest)@map %>% leaflet::clearBounds()
tforest_points <- mapedit::editMap(map_tforest)
tf_db <- tforest_points$drawn['geometry']
tf_db$type = "Tropical Forest"
cloudsen2 <- rbind(cloudsen2, tf_db)



## Temparated Forest
tdforest <- landuse_world$eq(3)
map_tdforest <- Map$addLayer(tdforest$updateMask(tdforest), list(min = 1, max = 10), name = "tdf") +
  map_base
map_tdforest <- mapview(cloudsen2, map_tdforest)@map %>% leaflet::clearBounds()
tdforest_points <- mapedit::editMap(map_tdforest)
tdf_db <- tdforest_points$drawn['geometry']
tdf_db$type = "Temparated Forest"
cloudsen2 <- rbind(cloudsen2, tdf_db)


## Grass/Crop
crop <- landuse_world$eq(4)
map_grass <- Map$addLayer(crop$updateMask(crop), list(min = 1, max = 10), name = "crop") +
  map_base
map_grass <- mapview(cloudsen2, map_grass)@map %>% leaflet::clearBounds()
grass_points <- mapedit::editMap(map_grass)
grass_db <- grass_points$drawn['geometry']
grass_db$type = "Grass/Crop"
cloudsen2 <- rbind(cloudsen2, grass_db)


## shrubland
shrubland <- landuse_world$eq(5)
map_shrubland <- Map$addLayer(shrubland$updateMask(shrubland), list(min = 1, max = 10), name = "shrubland") +
  map_base
map_shrubland <- mapview(cloudsen2, map_shrubland)@map %>% leaflet::clearBounds()
shrubland_points <- mapedit::editMap(map_shrubland)
shrubland_db <- shrubland_points$drawn['geometry']
shrubland_db$type = "Shrubland"
cloudsen2 <- rbind(cloudsen2, shrubland_db)


## Snow
snow <- landuse_world$eq(6)
map_snow <- Map$addLayer(snow$updateMask(snow), list(min = 1, max = 10), name = "snow") +
  map_base
map_snow <- mapview(cloudsen2, map_snow)@map %>% leaflet::clearBounds()
snow_points <- mapedit::editMap(map_snow)
snow_db <- snow_points$drawn['geometry']
snow_db$type = "Snow"
cloudsen2 <- rbind(cloudsen2, snow_db)

## Urban
urban <- landuse_world$eq(7)
map_urban <- Map$addLayer(urban$updateMask(urban), list(min = 1, max = 10), name = "urban") +
  map_base
map_urban <- mapview(cloudsen2, map_urban)@map %>% leaflet::clearBounds()

urban_points <- mapedit::editMap(map_urban)
urban_db <- urban_points$drawn['geometry']
urban_db$type = "Urban"
cloudsen2 <- rbind(cloudsen2, urban_db)


## Water
water <- landuse_world$eq(8)
map_water <- Map$addLayer(water$updateMask(water), list(min = 1, max = 10), name = "water") +
  map_base
map_water <- mapview(cloudsen2, map_water)@map %>% leaflet::clearBounds()

water_points <- mapedit::editMap(map_water)
water_db <- water_points$drawn['geometry']
water_db$type = "Water"
cloudsen2 <- rbind(cloudsen2, water_db)


## Wetlands
wetlands <- landuse_world$eq(9)
map_wetlands <- Map$addLayer(wetlands$updateMask(wetlands), list(min = 1, max = 10), name = "wetlands") +
  map_base
map_wetlands <- mapview(cloudsen2, map_wetlands)@map %>% leaflet::clearBounds()
wetlands_points <- mapedit::editMap(map_wetlands)
wetlands_db <- wetlands_points$drawn['geometry']
wetlands_db$type = "Wetlands"
cloudsen2 <- rbind(cloudsen2, wetlands_db)

write_sf(cloudsen2, "data/cloudsen2.geojson")
