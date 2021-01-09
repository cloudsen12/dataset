library(googlesheets4)
library(googledrive)
library(tidyverse)
library(jsonlite)
library(rgee)
library(sf)
source("src/utils.R")

# 1. Initialize Earth Engine
ee_Initialize("csaybar", drive = TRUE, gcs = TRUE)


# 2. List all metadata files
jsonfiles <- drive_metadata_json()
metadata <- "metadata_0019.json"

# 3. Generate the db
update_cloudsen12 <- function(jsonfiles, metadata) {
  # 3.1 point id
  point_id <- metadata
  point_number <- gsub("[a-z]|_|\\.","",point_id) %>% as.numeric()

  # 3.2 Assign point to a  labeler
  labeler <- assign_imgs(point_number)[point_number]

  # 3.3 Read metadata JSON
  jsonfile_r <- search_metajson(pattern = metadata, clean = FALSE) %>% read_json()

  # 3.4 SENTINEL-2 ID
  s2_ids <- names(jsonfile_r)[1:5]

  # 3.5 SENTINEL-1 ID
  s1_ids <- names(jsonfile_r)[1:5]

  # 3.6 GET Point
  # st_point <- st_sfc(geometry = st_point(c(jsonfile_r$x, jsonfile_r$y)), crs = 4326)
  # crs_kernel <- ee$Image(sprintf("COPERNICUS/S2/%s", s2_ids[1]))$select(0)$projection()$getInfo()$crs




}



# Sentinel2 ID

# POLYGON


save(drive_jsonfile, file = drive)

load(drive)
drive_jsonfile

# read metadata

point_n <- nrow(drive_jsonfile)


# 3. Generate area



dataset <- read_sheet("https://docs.google.com/spreadsheets/d/1LpW9JY2BdhlQvAObD1BCzoBiliNRnU3fCRWMQJFRvoM/edit#gid=0")





new_points <- c(39, 610,451,221,654)
for (new_point in new_points) {

  dataset_creator_chips(
    jsonfile = jsonfile,
    output_final = "/home/csaybar/Desktop/cloudsen12"
  )
}



dataset$id <- 10

sheet_write(
  data = dataset,
  ss = "https://docs.google.com/spreadsheets/d/1LpW9JY2BdhlQvAObD1BCzoBiliNRnU3fCRWMQJFRvoM/edit#gid=0",
  sheet = "main"
)
