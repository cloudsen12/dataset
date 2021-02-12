#' Script to create all the images in the cloudsen12 dataset
#' @author csaybar

# 1. Libraries
library(googleCloudStorageR)
library(googledrive)
library(tidyverse)
library(jsonlite)
library(mapview)
library(mapedit)
library(raster)
library(scales)
library(stars)
library(purrr)
library(grid)
library(rgee)
library(png)
library(sf)
library(sp)
source("src/utils.R")

# 2. Initialize Earth Engine
ee_Initialize("aybar1994", drive = TRUE, gcs = TRUE)

# 3. Load points with desired cloud average (after run point_creator.R)
local_cloudsen2_points <- read_sf("data/cloudsen2_potential_points.geojson")

# 4. Classify (label) images in clear, almost clear, low-cloudy, mid-cloudy, cloudy
# for (index in 1451:nrow(local_cloudsen2_points)) {
#   cloudsen2_row <- local_cloudsen2_points[index,]
#   select_dataset_thumbnail_creator(cloudsen2_row = cloudsen2_row)
# }

# 5. List all the metadata
jsonfile <- search_metajson(pattern = "metadata_0019.json", clean = FALSE)


# 6. Download all images in IRIS format :)
dataset_creator_chips(
  jsonfile = jsonfile,
  #upgrade_db = FALSE,
  output_final = "/home/csaybar/Desktop/cloudsen12"
)


# 7. Calibration
# metadata_f <- list.files("metadata/",pattern = "\\.json$",full.names = TRUE)
# detected_points_t <- sapply(metadata_f, detect_points)
# detected_points <- names(detected_points_t[detected_points_t >= 1])
# for (detected_point in detected_points) {
#   dataset_creator_chips(
#     jsonfile = detected_point,
#     #upgrade_db = FALSE,
#     output_final = "/home/csaybar/Desktop/cloudsen12"
#   )
# }


