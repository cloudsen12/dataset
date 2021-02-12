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
ee_Initialize()

local_cloudsen2_points <- read_sf("data/cloudsen2_potential_points.geojson")
extra_LC <- cloudsen12_lc()


landuse_01 <- ee_extract(extra_LC, local_cloudsen2_points$geometry[1:5000], scale = 100)
landuse_02 <- ee_extract(extra_LC, local_cloudsen2_points$geometry[5001:10000], scale = 100)
landuse_03 <- ee_extract(extra_LC, local_cloudsen2_points$geometry[10001:12023], scale = 100)
landuse_total <- c(landuse_01$land_cover, landuse_02$land_cover, landuse_03$land_cover)

table_sum <- na.omit(local_cloudsen2_points[c("type", "value")]) %>% st_drop_geometry() %>% unique() %>% arrange(value)
local_cloudsen2_points$value <- landuse_total

#landuse
landuse_factor <- as.factor(local_cloudsen2_points$value)
levels(landuse_factor) <- table_sum$type
local_cloudsen2_points$type <- landuse_factor

# pallete
pallete <- (ee$ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global") %>% ee_get(4))$first()$
  get("discrete_classification_class_palette")$
  getInfo()
landuse_factor <- as.factor(local_cloudsen2_points$value)
levels(landuse_factor) <- pallete
local_cloudsen2_points$palette <- landuse_factor


# probabilities
local_cloudsen2_points <- local_cloudsen2_points %>% arrange(id)
set.seed(100)
new_probs <- gen_rcloudpoints(length(5851:12023))
local_cloudsen2_points[5851:12023,4:8] <- new_probs

local_cloudsen2_points$id[5851:12023] <- 5851:12023
write_sf(local_cloudsen2_points, "data/cloudsen2_potential_points.geojson", delete_dsn = TRUE)


# mapviewOptions(fgb = FALSE)
# m1 <- mapview(local_cloudsen2_points, zcol = "type")
# mapshot(m1, "/home/csaybar/Documents/Github/cloudsen12/map/index.html")

#---------------------------------------------------------------------------

# 2. Initialize Earth Engine
ee_Initialize("aybar1994", drive = TRUE, gcs = TRUE)

# 3. Load points with desired cloud average (after run point_creator.R)
local_cloudsen2_points <- read_sf("data/cloudsen2_potential_points.geojson")


# Read metadata final
drive_jsonfile <- drive_ls(path = as_id("1fBGAjZkjPEpPr0p7c-LtJmfbLq3s87RK"))
full_points <- sprintf("metadata_%04d.json",1:1449)
# gsub("metadata_|\\.json$","", drive_jsonfile$name) %>% as.numeric() %>% max()

# upgrade good/bad field
bad_points <- full_points[!full_points %in% drive_jsonfile$name] %>%
  gsub("metadata_|\\.json$","",.) %>%
  as.numeric()
local_cloudsen2_points$good[1:1449] <- TRUE
local_cloudsen2_points$good[bad_points] <- FALSE

# Randomly split cloudsen12ms
# high-quality labeling -> 2400 (- 1449) = 951
# scribble labeling -> 4800
# no annotation -> 4800
set.seed(100)
select_points <- 1450:12023
hq_labels <- sample(select_points, 974)
set.seed(100)
s_labels <- sample(select_points[!select_points %in% hq_labels], 4800)
set.seed(100)
n_labels <- select_points[!select_points %in% c(hq_labels,s_labels)]


# point hq annotation
hq_points <- local_cloudsen2_points[c(1:1449, hq_labels), ]
hq_points$label <- "high_quality"
hq_points$labeler <- sprintf("labeler_%02d", c(rep(1:3, nrow(hq_points)/3), 1, 2))


# scribble annotation
s_points <- local_cloudsen2_points[s_labels, ]
s_points$label <- "scribble_annotation"
s_points$labeler <- sprintf("labeler_%02d", rep(1:3, nrow(s_points)/3))

# no annotation
n_points <- local_cloudsen2_points[n_labels, ]
n_points$label <- "no_annotation"
n_points$labeler <- "no_label"

final_dataset_points <- rbind(hq_points, s_points, n_points)
write_sf(final_dataset_points, "data/cloudsen2_potential_points.geojson", delete_dsn = TRUE)

#. ..............................
library(rgee)
library(sf)

sample_labeler <- function(n) {
  if (n == 1) {
    n + sample(1:2, 1)
  } else if(n == 2) {
    n + sample(c(-1, 1), 1)
  } else if(n == 3) {
    n + sample(c(-1, -2), 1)
  }  else if(n == 4) {
    n
  } else {
    stop("n")
  }
}

ee_Initialize()

cloudsen12 <- read_sf("data/cloudsen2_potential_points.geojson")
bynames <- cloudsen12$labeler %>% as.factor()
levels(bynames) <- 1:4
validator <- sapply(as.numeric(bynames), sample_labeler) %>% as.factor()
levels(validator) <- c("Jhomira", "Eduardo", "Fernando", "no_label")
levels(bynames) <- c("Jhomira", "Eduardo", "Fernando", "no_label")
cloudsen12$labeler <- bynames %>% as.character()
cloudsen12$validator <- validator %>% as.character()

write_sf(cloudsen12, "data/cloudsen2_potential_points.geojson",delete_dsn = TRUE)
