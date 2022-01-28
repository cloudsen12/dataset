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

ee_Initialize()

# 2. Python Libraries
ee_cloud <- import("ee_ipl_uv")
rio <- import("rasterio")
np <- import("numpy")

setwd("/home/csaybar/Documents/Github/cloudsen12/dataset/src/")
source_python("main.py")

ee_Initialize("csaybar", drive = TRUE, gcs = TRUE)

# 3. Auxiliary functions
setwd("/home/csaybar/Documents/Github/cloudsen12/dataset/")
source("src/utils.R")

# Identify points
googlesheets4::gs4_auth("s1078735@stud.sbg.ac.at")


# -----------------------------------------------------------------
dataset_creator_chips <- function(s2_id,
                                  point,
                                  kernel_size = c(255, 255),
                                  output_final = "cloudsen12/") {

  # 1. Read JSON file
  s2_id <- sprintf("COPERNICUS/S2/%s", s2_id)

  # 2. Create a st_point representing the center of the tile (255x255)
  bpoints <- list.files("/media/csaybar/Elements/cloudSEN12/", full.names = TRUE)
  allpoints <- lapply(bpoints, function(x) list.files(x, full.names = TRUE)) %>% unlist()
  detect_point <- allpoints[basename(allpoints) %in%  tolower(basename(point))]
  raster_ref <- list.files(detect_point, "sen2cor", full.names = TRUE, recursive = TRUE)[1] %>%
    raster()
  st_point <- st_sfc(st_point(c(mean(xFromCol(raster_ref)), mean(yFromRow(raster_ref)))), crs = crs(raster_ref)@projargs)
  crs_kernel <- sprintf("EPSG:%s", st_crs(raster_ref)$epsg)
  ee_point <- ee$Geometry$Point(st_point[[1]], proj = crs_kernel)

  # 4.1 S2 ID and dates
  s2_img <- ee$Image(s2_id)
  s2_date <- ee_get_date_img(s2_img)[["time_start"]]

  # 4.2 S1 ID
  s1_id <- ee_get_s1(point = ee_point, s2_date = s2_date)
  s1_img <- ee$Image(s1_id)
  s1_date <- ee_get_date_img(s1_img)[["time_start"]]

  # 4.3 Create an Image collection with S2, S1 and cloud mask information
  s2_s1_img <- ee_merge_s2_full(s2_id, s1_id, s2_date)

  # 4.4 Add the shadow direction
  s2_fullinfo <- s2_s1_img %>%
    ee$Image$addBands(shadow_direction(image = s2_img))

  # 4.5 Create a 511x511 tile (list -> data_frame -> sp -> raster)
  band_names <- c(s2_fullinfo$bandNames()$getInfo(), "x", "y")
  s2_img_array <- s2_fullinfo$addBands(s1_img) %>%
    ee$Image$addBands(ee$Image$pixelCoordinates(projection = crs_kernel)) %>%
      ee$Image$neighborhoodToArray(
        kernel = ee$Kernel$rectangle(kernel_size[1], kernel_size[2], "pixels")
      ) %>%
      ee$Image$sampleRegions(ee$FeatureCollection(ee_point),
                             projection = crs_kernel,
                             scale = 10) %>%
      ee$FeatureCollection$getInfo()

  if (length(s2_img_array$features) == 0) {
    stop("The datacube doesn't have values ... Probably a bug in S1")
  }

  extract_fn <- function(x) as.numeric(unlist(s2_img_array$features[[1]]$properties[x]))
  image_as_df <- do.call(cbind,lapply(band_names, extract_fn))
  colnames(image_as_df) <- band_names
  image_as_tibble <- as_tibble(image_as_df)
  coordinates(image_as_tibble) <- ~x+y
  sf_to_stack <- function(x) rasterFromXYZ(image_as_tibble[x])
  final_stack <- stack(lapply(names(image_as_tibble), sf_to_stack))
  crs(final_stack) <- st_crs(crs_kernel)$proj4string

  # 4.8 Add 360 if shadow direction is negative
  final_stack[[24]] <- fix_shadow_direction(final_stack[[24]])

  # 4.9 Prepare folders for iris
  output_final <- "/home/csaybar/Desktop/fixing/initial/"
  output_final_d <- sprintf("%s/dataset", output_final)
  point_name <- tolower(point)
  output_final_folder <- sprintf("%s/dataset/%s/%s", output_final, point_name, basename(s2_id))

  metadata_main <- sprintf("%s/dataset/%s/cloud_segmentation_%s.json", output_final, point_name, point_name)
  metadata_spec <- sprintf("%s/dataset/%s/%s/metadata.json", output_final, point_name, basename(s2_id))

  dir.create(sprintf("%s/input", output_final_folder), showWarnings = FALSE, recursive = TRUE)
  dir.create(sprintf("%s/target", output_final_folder), showWarnings = FALSE, recursive = TRUE)
  dir.create(sprintf("%s/thumbnails", output_final_folder), showWarnings = FALSE, recursive = TRUE)

  # 4.9 Create cloud-segmentation.json (main file in iris software)
  metadata_main <- sprintf("%s/dataset/%s/cloud_segmentation_%s.json", output_final, point_name, point_name)
  ee_create_cloudseg(path = metadata_main)

  # 4.11 Generate a script to upload results to Google Drive database
  generate_script(dirname(metadata_main))

  # 4.12 Save all features inside the input folder
  bandnames <- c(paste0("B",1:8), "B8A", paste0("B", 9:12), "CDI", "VV", "VH", "angle", "elevation", "landuse", "cloudshadow_direction")

  # 4.14 Save input values
  input_data <- raster::stack(
    final_stack[[1:13]]/10000, final_stack[[14]], final_stack[[15:17]], final_stack[[22:24]]
  )
  input_spec <- sprintf("%s/input/%s.tif", output_final_folder, bandnames)
  lapply(1:20, function(x) writeRaster(input_data[[x]], input_spec[x], overwrite = TRUE))

  # B10 threshold
  input_spec_b10 <- sprintf("%s/input/%s.tif", output_final_folder, "B10_threshold")
  writeRaster((input_data[[11]] > 0.03), input_spec_b10, overwrite = TRUE)

  if (maxValue(final_stack[[14]] == -99)) {
    stop("Sentinel2 with NaN data")
  }

  # 4.15 Save target values
  # 18-19 -> cmask_s2cloudness| cmask_s2cloudness_reclass (0,1)
  # 20-21 -> cmask_sen2cor | cmask_sen2cor_reclass (0,1,2)
  # 25 -> IPL_cloudmask_reclass
  create_target_raster(
    final_stack = final_stack,
    output_final_folder = output_final_folder
  )

  # 4.17 Create a thumbnail
  mapview::mapview(final_stack)
  ee_generate_thumbnail(
    s2_id = s2_id,
    final_stack =  final_stack,
    crs_kernel =  crs_kernel,
    output_final_folder = output_final_folder
  )
}


# Functions
cloudsen12_up_local <- function(point_dir, fmask_files) {
  message("Processing: ", basename(point_dir))

  # 1. Identify the point :)
  point_dir_full_files <- list.files(point_dir,full.names = TRUE)

  # 2. get only the folder with s2/s1 data :D
  files_points_only_img_folder <- point_dir_full_files[grepl("^[0-9]", basename(point_dir_full_files))]
  if (length(files_points_only_img_folder) != 1) {
    stop("Initial regex condition fail :/")
  }

  # 3. Load raster ref
  refraster <- list.files(
    path = files_points_only_img_folder[1],
    pattern = "sen2cor_real",
    recursive = TRUE,
    full.names = TRUE
  ) %>% raster()

  # 4. QA60 generation
  for (index in 1:1) {
    fullname_s2f <- files_points_only_img_folder[index]
    q60 <- QA60_creator(basename(fullname_s2f), refraster)
    q60_out <- sprintf("%s/target/qa60.tif", fullname_s2f)
    writeRaster(x = q60, filename = q60_out, overwrite = TRUE)
  }

  # 5. Fmask generation
  for (index in 1:1) {
    fullname_s2f <- files_points_only_img_folder[index]
    rfmask <- FMask_creator(basename(fullname_s2f), refraster, fmask_files)
    rfmask_out <- sprintf("%s/target/fmask.tif", fullname_s2f)
    writeRaster(x = rfmask, filename = rfmask_out, overwrite = TRUE)
  }
  # 6. UNET
}

QA60_creator <- function(s2_id, raster_ref) {
  # Get centroid coordinates
  cente_sp <- st_sfc(st_point(apply(coordinates(raster_ref), 2, mean)), crs = crs(raster_ref)@projargs)
  crs_sp <- st_crs(read_stars(raster_ref@file@name, proxy = TRUE))$epsg

  # Sentinel-2 Level 1C
  s2_img <- ee$Image(sprintf("COPERNICUS/S2/%s", s2_id))$select("QA60")

  # Create a st_point representing the center of the tile (255x255)
  crs_kernel <- sprintf("EPSG:%s", crs_sp)
  point_utm <- st_transform(cente_sp, crs_kernel)
  ee_point <- ee$Geometry$Point(point_utm[[1]], proj = crs_kernel)
  s2_img <- ee$Image$reproject(s2_img, crs_kernel)

  # 4. Create a 511x511 tile (list -> data_frame -> sp -> raster)
  s2_img %>%
    ee$Image$addBands(ee$Image$pixelCoordinates(projection = crs_kernel)) %>%
    ee$Image$neighborhoodToArray(
      kernel = ee$Kernel$rectangle(255, 255, "pixels")
    ) %>%
    ee$Image$sampleRegions(ee$FeatureCollection(ee_point),
                           projection = crs_kernel,
                           scale = 10) %>%
    ee$FeatureCollection$getInfo() -> s2_img_array

  band_names <- c("QA60", "x", "y")
  extract_fn <- function(x) as.numeric(unlist(s2_img_array$features[[1]]$properties[x]))
  image_as_df <- do.call(cbind,lapply(band_names, extract_fn))
  colnames(image_as_df) <- band_names
  image_as_tibble <- as_tibble(image_as_df)
  coordinates(image_as_tibble) <- ~x+y
  sf_to_stack <- function(x) rasterFromXYZ(image_as_tibble[x])
  final_stack <- stack(lapply(names(image_as_tibble), sf_to_stack))
  crs(final_stack) <- st_crs(crs_kernel)$proj4string
  final_stack
}

FMask_creator <- function(s2_id, raster_ref, fmask_files, tempfile = "/home/csaybar/z_.tif") {
  allfmaskfiles_b <- basename(fmask_files)
  id <- ee$Image(sprintf("COPERNICUS/S2/%s", s2_id))$get("GRANULE_ID")$getInfo()
  fmask_file <- fmask_files[allfmaskfiles_b %in% sprintf("%s_Fmask4.tif", id)][1]

  crs <- st_crs(raster_ref)$epsg
  rextent <- extent(raster_ref)
  xmin <- rextent[1]
  xmax <- rextent[2]
  ymin <- rextent[3]
  ymax <- rextent[4]
  system(
    sprintf(
      "gdalwarp '%s' '%s' -overwrite -te %s %s %s %s -q -tr 10 10 -t_srs %s",
      fmask_file, tempfile,
      xmin, ymin,
      xmax, ymax, sprintf("EPSG:%s", crs)
    )
  )
  raster(tempfile)
}


Rreshape <- function(x, output, raster_ref) {
  crs <- st_crs(raster_ref)$epsg
  rextent <- extent(raster_ref)
  xmin <- rextent[1]
  xmax <- rextent[2]
  ymin <- rextent[3]
  ymax <- rextent[4]
  system(
    sprintf(
      "gdalwarp '%s' '%s' -overwrite -co 'COMPRESS=LZW' -te %s %s %s %s -q -tr 10 10 -t_srs %s",
      x, output,
      xmin, ymin,
      xmax, ymax, sprintf("EPSG:%s", crs)
    )
  )
  raster(output)
}


remove_borders <- function(x) {
  xmin <- extent(x)[1] + 10
  xmax <- extent(x)[2] - 10
  ymin <- extent(x)[3] + 10
  ymax <- extent(x)[4] - 10
  rr <- raster(
    nrows = 509, ncols = 509,
    xmn = xmin, xmx = xmax,
    ymn = ymin, ymx = ymax,
    crs = crs(x)
  )
  rr[] <- x[2:510,2:510]
  rr
}

create_npy_file <- function(fpath, npy_output) {
  files_input_ffolder <- list.files(fpath, "\\.tif$",full.names = TRUE)
  files_input_ffolder_nob10 <- files_input_ffolder[!basename(files_input_ffolder) %in% "B10_threshold.tif"]
  names(files_input_ffolder_nob10) <- gsub("\\.tif$", "", basename(files_input_ffolder_nob10))

  if (length(files_input_ffolder_nob10) != 20) {
    stop("create_npy_file error")
  }

  # Names order
  names_in_order <- c(
    "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11",
    "B12", "VH", "VV", "angle", "CDI", "cloudshadow_direction", "elevation",
    "landuse"
  )

  Rmatrix <- files_input_ffolder_nob10[names_in_order] %>%
    read_stars() %>%
    merge() %>%
    '[['("X") %>%
    '['(2:510, 2:510, ) %>%
    '*'(10000) %>%
    np$moveaxis(2L, 0L) %>%
    py$np_storage(npy_output)
  npy_output
}



thumbnails_creator <- function(thumbnail_in, thumbnail_out, raster_ref) {
  # Get the CRS
  crs <- st_crs(raster_ref)$epsg

  # getting te and tr parameters
  r <- raster(raster_ref)
  centroid <- st_sfc(st_point(c(mean(xFromCol(r)), mean(yFromRow(r)))), crs = crs(r)@projargs)
  ref_b <- st_buffer(centroid, 475*10, endCapStyle="SQUARE")
  xmin <- as.numeric(st_bbox(ref_b))[1] - 5
  ymin <- as.numeric(st_bbox(ref_b))[2] - 5
  xmax <- as.numeric(st_bbox(ref_b))[3] + 5
  ymax <- as.numeric(st_bbox(ref_b))[4] + 5

  # save the resutls
  system(
    sprintf(
      "gdalwarp '%s' '%s' -overwrite -ot 'Byte'  -co 'COMPRESS=LZW' -te %s %s %s %s -tr 10 10 -q -t_srs %s",
      thumbnail_in, thumbnail_out,
      xmin, ymin, xmax, ymax, sprintf("EPSG:%s", crs)
    )
  )
}

new_db_migration_local <- function(point, output) {
  message("Processing point: ", basename(point))

  # 1. Identify the point :)
  point_dir_full_files <- list.files(point, full.names = TRUE)

  # 2. get only the folder with s2/s1 data :D
  files_points_only_img_folder <- point_dir_full_files[grepl("^[0-9]", basename(point_dir_full_files))]
  if (length(files_points_only_img_folder) != 1) {
    stop("Initial regex condition fail :/")
  }

  # 3. Create directories
  s2ids <- basename(files_points_only_img_folder)
  fpaths <- sprintf("%s/%s/%s/models", output, basename(point), s2ids)
  lapply(
    X = 1:1,
    FUN = function(x) dir.create(path = fpaths[x],
                                 recursive = TRUE,
                                 showWarnings = FALSE)
  )

  # 4. Create raster ref
  s2corfiles <- list.files(
    path = point,
    pattern = "sen2cor_real\\.tif$",
    full.names = TRUE,
    recursive = TRUE
  )
  s2cor_stack <- stack(s2corfiles) # does it create the stack?
  raster_ref <- s2cor_stack[[1]] %>% remove_borders


  # 5. Sen2COR
  for (index in 1:1) {
    s2xid <- basename(dirname(dirname(s2corfiles[index])))
    outf <- sprintf("%s/%s/%s/models/sen2cor.tif", output, basename(point), s2xid)
    Rreshape(s2corfiles[index], outf, raster_ref)
  }

  # 6. Sencloudness
  s2cloudlessfiles <- list.files(
    path = point,
    pattern = "s2cloudness_prob\\.tif$",
    full.names = TRUE,
    recursive = TRUE
  )
  for (index in 1:1) {
    s2xid <- basename(dirname(dirname(s2cloudlessfiles[index])))
    outf <- sprintf("%s/%s/%s/models/s2cloudless.tif", output, basename(point), s2xid)
    Rreshape(s2cloudlessfiles[index], outf, raster_ref)
  }

  # 7. qa60
  qafiles <- list.files(
    path = point,
    pattern = "qa60\\.tif$",
    full.names = TRUE,
    recursive = TRUE
  )
  for (index in 1:1) {
    s2xid <- basename(dirname(dirname(qafiles[index])))
    outf <- sprintf("%s/%s/%s/models/qa60.tif", output, basename(point), s2xid)
    Rreshape(qafiles[index], outf, raster_ref)
  }

  # 8. fmask
  fmaskfiles <- list.files(
    path = point,
    pattern = "fmask\\.tif$",
    full.names = TRUE,
    recursive = TRUE
  )
  for (index in 1:1) {
    s2xid <- basename(dirname(dirname(fmaskfiles[index])))
    outf <- sprintf("%s/%s/%s/models/fmask.tif", output, basename(point), s2xid)
    Rreshape(fmaskfiles[index], outf, raster_ref)
  }

  # 9. Create npy
  npoutfs <- rep(NA, 1)
  for (index in 1:1) {
    s2folder <- dirname(fpaths[index])
    s2npyid <- basename(s2folder)
    inf <-sprintf("%s/%s/input", point, s2npyid)
    outf <- sprintf("%s/input.npy", s2folder)
    create_npy_file(inf, outf)
    npoutfs[index] <- outf
  }


  # 10. Thumbnail
  for (index in 1:1) {
    s2folder <- dirname(fpaths[index])
    s2npyid <- basename(s2folder)
    thumbnail_in <-sprintf("%s/%s/thumbnails/thumbnail.tif", point, s2npyid)
    thumbnail_out <- sprintf("%s/thumbnail.tif", s2folder)
    thumbnails_creator(thumbnail_in, thumbnail_out, raster_ref)
  }

  # 11. Create UNET
  for (index in 1:1) {
    npy_file_input <- npoutfs[index]
    model_results01 <- py$unet_main(npy_file = npy_file_input, namemodel = "rgbiswir")
    model_results02 <- py$unet_main(npy_file = npy_file_input, namemodel = "rgbi")
    model_results_raster_01 <- raster_ref
    model_results_raster_01[] <- model_results01

    model_results_raster_02 <- raster_ref
    model_results_raster_02[] <- model_results02

    writeRaster(
      x = model_results_raster_01*100,
      filename = sprintf("%s/models/UNET_rgbiswir.tif", dirname(npy_file_input)),
      datatype = "INT1U",
      overwrite = TRUE
    )
    writeRaster(
      x = model_results_raster_02*100,
      filename = sprintf("%s/models/UNET_rgbi.tif", dirname(npy_file_input)),
      datatype = "INT1U",
      overwrite = TRUE
    )
  }

  # 12. Check
  all_files <- list.files(sprintf("%s/%s", output, basename(point)), full.names = TRUE,recursive = TRUE)
  thumnbnail_stk <- all_files[grep("thumbnail", all_files)]
  read_stars(thumnbnail_stk)
  nothumnbnail_stk <- all_files[!grepl("thumbnail|input\\.npy", all_files)]
  read_stars(nothumnbnail_stk)
}


# ------------------------------------------------------------------------
googlesheets4::gs4_auth("s1078735@stud.sbg.ac.at")
bad_good_db <- googlesheets4::read_sheet("1_5G5SBvdODn4PDKP_a-9RYhhtdVpQHyDKB8SWggP6-s", 13)


#check ids
xxx <- sapply(bad_good_db$NEW_ID, function(x) {
  tryCatch(
    expr = {
      ee$Image(sprintf("COPERNICUS/S2_SR/%s", x))$getInfo()
      TRUE
    },
    error = function(e) FALSE
  )}
)

FIXING_FOLDER <- "/home/csaybar/Desktop/fixing"
FMASK_FOLDER <- "/home/csaybar/Downloads/FMASK4/"


FMASK_FOLDER_FILES <- list.files(FMASK_FOLDER, "\\.tif$", recursive = TRUE, full.names = TRUE)
index <- 14
which(xxx)[24]
for (index in which(xxx)[24:168]) {
  point <- bad_good_db$Point[index] %>% tolower()
  s2_id <- bad_good_db$NEW_ID[index]
  dataset_creator_chips(s2_id, point, output_final = FIXING_FOLDER)
  point_dir <- sprintf("%s/initial/dataset/%s", FIXING_FOLDER, point)
  cloudsen12_up_local(point_dir = point_dir, fmask_files = FMASK_FOLDER_FILES)
  # -----------------------------------
  lpoint <- sprintf("%s/initial/dataset/%s", FIXING_FOLDER, basename(point))
  new_db_migration_local(point = lpoint, output = "/home/csaybar/Desktop/fixing/final/")
}


