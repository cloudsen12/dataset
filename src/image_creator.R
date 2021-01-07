#' Create image in cloudsen12
#' @author csaybar
#'
#' Script used to manually select images in cloudsen12.
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

set.seed(101)
source("src/utils.R")

ee_Initialize("aybar1994", drive = TRUE, gcs = TRUE)

# 1. Load points with desired cloud average (after run point_creator.R)
local_cloudsen2_points <- read_sf("data/cloudsen2_potential_points.geojson")

# 2. Classify images in clear, almost clear, low-cloudy, mid-cloudy, cloudy
for (index in 1451:nrow(local_cloudsen2_points)) {
  cloudsen2_row <- local_cloudsen2_points[index,]
  select_dataset_thumbnail_creator(
    cloudsen2_row = cloudsen2_row,
    n_images = "max",
    kernel_size = c(255, 255),
    data_range = c("2018-01-01", "2020-07-31"),
    output = "results/"
  )
}

# 3. Download images!
jsonfiles <- list.files(
  path = "/home/csaybar/Documents/Github/cloudsen12/dataset/results/",
  pattern = "\\.json$",
  recursive = TRUE
)

for (index in seq_alon(jsonfiles)) {
  dataset_creator_chips(
    jsonfile = jsonfile,
    kernel_size = c(255, 255),
    output = "results/"
  )
}
