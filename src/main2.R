#' Add QA60 and FMASK to cloudSEN12 folder
library(reticulate)
library(raster)
library(stars)
library(dplyr)
library(rgee)
library(sf)

ee_Initialize("aybar1994", drive = TRUE, gcs = TRUE)

np <- import("numpy")
source_python("main.py")


# Global parameters
FMASK_FOLDER <- "/home/csaybar/Downloads/FMASK4/"
FMASK_FOLDER_FILES <- list.files(FMASK_FOLDER, "\\.tif$", recursive = TRUE, full.names = TRUE)

# RUN 

## high
point_dirs <- list.files("/media/csaybar/Elements SE/cloudSEN12/high/", full.names = TRUE)
for (index in 1:2000) {
  print(index)
  cloudsen12_up_local(point_dirs[index], FMASK_FOLDER_FILES)  
}

## scribble
point_dirs <- list.files("/media/csaybar/Elements SE/cloudSEN12/scribble/", full.names = TRUE)
for (index in 1:2000) {
  print(index)
  cloudsen12_up_local(point_dirs[index], FMASK_FOLDER_FILES)  
}

## nolabel
point_dirs <- list.files("/media/csaybar/Elements SE/cloudSEN12/nolabel/", full.names = TRUE)
for (index in 1:6000) {
  print(index)
  cloudsen12_up_local(point_dirs[index], FMASK_FOLDER_FILES)  
}



img <- ee$Image("COPERNICUS/S2/20200625T153619_20200625T153621_T17MNP")
Map$centerObject(img)

mm <- Map$addLayer(img, list(min=0,max=4000))+
  Map$addLayer(img[["QA60"]]>0, list(min=0,max=1))

library(mapview)
mapview(rrr, mm)
rrr <- raster("/media/csaybar/Elements SE/cloudSEN12/high/point_0005/20200625T153619_20200625T153621_T17MNP/target/qa60.tif")
FMASK_FOLDER_FILES

id <- img$get("GRANULE_ID")$getInfo()
fmask_file <- fmask_files[allfmaskfiles_b %in% sprintf("%s_Fmask4.tif", id)][1]
file.copy(
  fmask_file,
  "/home/csaybar/Desktop/ff.tif"
)
