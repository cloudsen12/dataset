#' Create image in cloudsen12
#' @author csaybar
#'
#' Script used to manually select images in cloudsen12.

library(tidyverse)
library(jsonlite)
library(mapview)
library(mapedit)
library(raster)
library(scales)
library(stars)
library(grid)
library(rgee)
library(png)
library(sf)
library(sp)

set.seed(101)
source("src/utils.R")
ee_Initialize("csaybar", drive = TRUE)

# 1. Load points with desired cloud average (after run point_creator.R)
local_cloudsen2_points <- read_sf("data/cloudsen2_points.geojson")


# 2. Classify images in clear, almost clear, low-cloudy, mid-cloudy, cloudy
for (index in seq_len(nrow(local_cloudsen2_points))) {
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
for (index in seq_len(nrow(local_cloudsen2_points))) {
  cloudsen2_row <- local_cloudsen2_points[index,]
  dataset_creator_chips2(
    cloudsen2_row = cloudsen2_row,
    kernel_size = c(255, 255),
    output = "results/"
  )
}
