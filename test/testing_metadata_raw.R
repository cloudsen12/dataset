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
ee_Initialize("valeria")

source("src/utils.R")
ee_cloud <- import("ee_ipl_uv")



test_metadata_raw <- function() {
  dir_p  <- paste0(tempdir(), "/cd26ed5dc626f11802a652e81d02762e_s1078735@stud.sbg.ac.at")
  #download privileges
  download.file(
    url = "https://drive.google.com/uc?id=1dg02Ue6rkFa7xn7TU57bWu_gQxuxpGdY&export=download",
    destfile = dir_p
  )
  
  drive_auth("s1078735@stud.sbg.ac.at", token = dir_p)
  
  # Points Google Drive ID
  files_points_general <- googledrive::drive_ls(
    path = as_id("1NOxjbbtiyz2UqiJAJz6IsCNaLxOJlErr")
  )
  metadata_raw <- sprintf("point_%04d", 1451:12000)
  metadata_rawf <- metadata_raw[which(!metadata_raw %in% sort(files_points_general$name))]
  metadata_rawf
}

metadata_rawf <- test_metadata_raw()
# 3. Load points with desired cloud average (after run point_creator.R)
local_cloudsen2_points <- read_sf("data/cloudsen2_potential_points.geojson")

# 4. Classify (label) images in clear, almost clear, low-cloudy, mid-cloudy, cloudy
# cesar <- 5851:7909
# roy <- 7910:9968
# prudencio <- 9969:12023
lost_points <- gsub("point_", "", metadata_rawf) %>% as.numeric
lost_points_download <- lost_points[100:130]

select_dataset_thumbnail_creator_batch <- function(points, local_cloudsen2_points) {
  for (index in points) {
    print(index)
    cloudsen2_row <- local_cloudsen2_points[index,]
    try(select_dataset_thumbnail_creator(cloudsen2_row = cloudsen2_row))
  }
}

select_dataset_thumbnail_creator_batch(lost_points_download, local_cloudsen2_points)

