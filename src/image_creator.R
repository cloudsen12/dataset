#' Script to create all the images in the cloudsen12 dataset
#'
#' human label:
#'   0 -> clear
#'   1 -> thick cloud
#'   2 -> thin cloud
#'   3 -> cloud shadow
#'
#' @author csaybar

# 1. Libraries
library(googleCloudStorageR)
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

# 2. Initialize Earth Engine
ee_Initialize("csaybar", drive = TRUE, gcs = TRUE)

source("src/utils.R")
ee_cloud <- import("ee_ipl_uv")


# 3. Load points with desired cloud average (after run point_creator.R)
local_cloudsen2_points <- read_sf("data/cloudsen2_potential_points.geojson")

# 4. Classify (label) images in clear, almost clear, low-cloudy, mid-cloudy, cloudy
# cesar <- 5851:7909
# roy <- 7910:9968
# prudencio <- 9969:12023
select_dataset_thumbnail_creator_batch(cesar, local_cloudsen2_points)

# # # 5. Download all images
download_all_images(1:1500, local_cloudsen2_points)
