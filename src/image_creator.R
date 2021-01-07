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
for (index in 1451:nrow(local_cloudsen2_points)) {
  cloudsen2_row <- local_cloudsen2_points[index,]
  select_dataset_thumbnail_creator(cloudsen2_row = cloudsen2_row)
}

# 5. List all the metadata
jsonfiles <- list.files(
  path = "/home/csaybar/Documents/Github/cloudsen12/dataset/results/",
  pattern = "\\.json$",
  recursive = TRUE
)

# 6. Download all images in IRIS format :)
for (jsonfile in seq_alon(jsonfiles)) {
  dataset_creator_chips(jsonfile = jsonfile)
}
