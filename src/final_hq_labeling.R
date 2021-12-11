library(png)
library(doMC)
library(raster)
library(ggplot2)
library(foreach)
library(lubridate)
library(tidyverse)
library(googledrive)
library(googlesheets4)

googledrive::drive_auth("s1078735@stud.sbg.ac.at")
googlesheets4::gs4_auth("s1078735@stud.sbg.ac.at")

DATASET_FOLDER <- "/media/csaybar/Elements/cloudSEN12x/high/"
DATASET_OUTPUT <- "/home/csaybar/Desktop/labels_results/"

registerDoMC(10)

# FUNCTIONS --------------------------------------------------------------------
png_to_tif <- function(name, alpha = 1) {
  png_mtx <- png::readPNG(name)
  png_mtx <- png_mtx[2:510, 2:510, 1:3]
  rgb_img <- rgb(png_mtx[,,1], png_mtx[,,2], png_mtx[,,3])
  rgb_img[rgb_img == "#FFFFFF"] = 0
  rgb_img[rgb_img == "#FFFF00"] = 1
  rgb_img[rgb_img == "#00FF00"] = 2
  rgb_img[rgb_img == "#FF0000"] = 3
  rgb_img <- as.numeric(rgb_img)
  dim(rgb_img) <- c(509, 509)
  point_n <- strsplit(basename(name), "__")[[1]][1]
  s2id_n <- strsplit(basename(name), "__")[[1]][2]
  ref_r <- raster(sprintf("%s/%s/%s/models/sen2cor.tif", DATASET_FOLDER, point_n, s2id_n))
  ref_r[] <- rgb_img
  ref_r
}

hq_dataset <- function(point) {
  require <- list.files(sprintf("%s/%s", DATASET_FOLDER, point))
  files_found_it <- allfiles[gsub(".*__(.*)__manual\\.png$","\\1", basename(allfiles)) %in% require]
  files_found_it <- files_found_it[grepl(point, files_found_it)]
  if (length(files_found_it) != 5) {
    message(paste0(files_found_it, collapse = "\n"))
    stop("no 5")
  }
  lapply(files_found_it, function(x) {
    rr <- png_to_tif(x)
    xname <- gsub("\\.png", "\\.tif", basename(x))
    writeRaster(rr, sprintf("%s/%s", DATASET_OUTPUT, xname), overwrite = TRUE)
  })
  invisible(TRUE)
}
# ------------------------------------------------------------------------------


# 1. Obtain all the zip_files
drive_jsonfile <- drive_ls(
  path = as_id("133GLIclZiNRqBfLA7aEInXyOI-sj1bYB"),
  recursive = TRUE
)

# 2. Download all the zip files
drive_jsonfile_02 <- drive_jsonfile[grepl("\\.zip$", drive_jsonfile$name),]
tmpzipfile <- tempfile(fileext = ".zip")
tmpzfolder <- tempfile()
dir.create(tmpzfolder)
for (index in 1:nrow(drive_jsonfile_02)) {
  drive_download(
    file = drive_jsonfile_02[index,],
    path = tmpzipfile, 
    overwrite = TRUE
  )
  unzip(tmpzipfile, exdir = tmpzfolder)
}


#### Create TIFF files
allfiles <- list.files(tmpzfolder, "\\.png$", recursive = TRUE, full.names = TRUE)
points <- list.files(DATASET_FOLDER)


foreach(index = 1:length(points)) %dopar% {
  x <- hq_dataset(point = points[index])
}


zip("/home/csaybar/Desktop/high_quality.zip", list.files(DATASET_OUTPUT, full.names = TRUE), flags = "-j")
