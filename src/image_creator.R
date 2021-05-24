#' Main script to construct cloudSEN12. This script support three-step:
#'
#' (1) select_dataset_thumbnail_creator_batch: Tile selection.
#' (2) download_images: Download all the points from cloudSEN12 (friendly with IRIS).
#' (3) db_migration: Create .npy and order cloud detection models
#' (4) stac_feature_creator: Create single-file .json with image metadata following the STAC Item specification
#' (5) stac_featurecollection_creator: Merge all item files (.json) in a single file (FeatureCollection).
#'
#' human label:
#'   0 -> clear
#'   1 -> thick cloud
#'   2 -> thin cloud
#'   3 -> cloud shadow
#'
#' @author csaybar

# 1. R Libraries
library(googleCloudStorageR)
library(googlesheets4)
library(filesstrings)
library(googledrive)
library(reticulate)
library(rasterVis)
library(tidyverse)
require(gridExtra)
library(jsonlite)
library(mapview)
library(mapedit)
library(raster)
library(scales)
library(mmand)
library(stars)
library(purrr)
library(grid)
library(Orcs)
library(rgee)
library(png)
library(sf)
library(sp)

# 2. Python Libraries
ee_cloud <- import("ee_ipl_uv")
rio <- import("rasterio")
np <- import("numpy")

# 3. Auxiliary functions
source("src/utils.R")

# 2. Initialize Earth Engine, GD and GS4 (the user must have permission to write in csaybar and roy GD folder)
# https://drive.google.com/drive/u/1/folders/1aTPIZ974zvtti6a02eiMyZaIf_Rp8QEc (ROY DRIVE)
# https://drive.google.com/drive/u/0/folders/1tcFgbP3SLovBy3UNs7OPMWNOO5PuiCuP (CSAYBAR DRIVE)
# https://console.cloud.google.com/storage/browser/cloudsen12 (CSAYBAR bucket)
ee_Initialize("csaybar", drive = TRUE, gcs = TRUE)
drive_auth("aybar1994@gmail.com")
gs4_auth("aybar1994@gmail.com")


# 3. Load Potential points (points created with point_creator.R)
local_cloudsen2_points <- read_sf("data/cloudsen2_potential_points.geojson")
# tail(local_cloudsen2_points$label[1:1450],10)

# 4. Tile selection on each point and image
select_dataset_thumbnail_creator_batch(
  points = 6707:6777,
  local_cloudsen2_points = local_cloudsen2_points,
  n_images = "max",
  kernel_size = c(255, 255),
  data_range = c("2018-01-01", "2020-07-31"),
  output = "/home/csaybar/cloudSEN12/"
)

# 5. Download all CLOUDSEN12 images with a format friendly with IRIS
lost_points <- gsub("metadata_|.json", "", basename(read.csv("/home/csaybar/Desktop/lost_metadata.csv")$metadata_lost)) %>% as.numeric()
download_cloudSEN12_images(
  points = lost_points,
  local_cloudsen2_points = local_cloudsen2_points,
  output = "/home/csaybar/Desktop/cloudsen12/"
)

# 6. Database migration, we restructured the database with a format easy to
#    ingest into a deep learning model.
to_migrate <- gsub("point_", "", list.files("/media/csaybar/Elements SE/cloudSEN12/high/")) %>% as.numeric()
fpoints <- list.files("/media/csaybar/Elements SE/cloudSEN12_f/high/", full.names = TRUE)
again_download <- list()
counter <- 1
for (index in seq_len(length(fpoints))) {
  length_p <- length(list.files(fpoints[index], recursive = TRUE))
  if (length_p != 25) {
    again_download[[counter]]  = index
    counter <- counter + 1
  }
}

# again_download %>% as.numeric()
db_migration_local_batch(
  points = to_migrate[again_download %>% as.numeric()],
  dataset_dir = "/media/csaybar/Elements SE/cloudSEN12/high/",
  output = "/media/csaybar/Elements SE/cloudSEN12_f/high/"
)

# db_migration_manual_labeling()
# db_migration_manual_fmask4()
# db_migration_manual_maja()
# db_migration_manual_unet()
# db_migration_manual_lightgbm()


# 7. Create STAC items (features) following the single file STAC Extension Specification
# https://github.com/stac-extensions/single-file-stac
stac_feature_creator_batch(
  points = 1,
  potencial_points = local_cloudsen2_points,
  output = "/home/csaybar/cloudSEN12/"
)

# 8. Change the container
local_container <- "/home/csaybar/cloudSEN12/points/"
gcs_container <- "https://storage.googleapis.com/cloudsen12/CloudSEN12/"
upgrade_link_gcs(json_file, gcs_container, local_container)



# 9. Merge all json in a same file (FeatureCollection)
json_list <- list.files("/home/csaybar/cloudSEN12/", "\\.json$", recursive = TRUE, full.names = TRUE)

# local_container <- "/home/csaybar/cloudSEN12/points/"
# gcs_container <- "https://storage.googleapis.com/cloudsen12/CloudSEN12/"
# upgrade_link_gcs(json_file, gcs_container, local_container)
# mapply(upgrade_link_gcs, json_list[1:5], MoreArgs = list(local_container=local_container, gcs_container=gcs_container))

stac_featurecollection_creator(
  json_list = json_list[1:5],
  outputfile = "/home/csaybar/Documents/Github/STAC/data/index.geojson"
)
