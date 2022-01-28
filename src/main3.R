#' Main script to construct cloudSEN12. This script support three-step:
library(reticulate)
library(raster)
library(stars)
library(dplyr)
library(stars)
library(rgee)
library(sf)

setwd("/home/csaybar/Documents/Github/cloudsen12/dataset/src/")
np <- import("numpy")
source_python("main.py")

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
  if (length(files_points_only_img_folder) != 5) {
    stop("Initial regex condition fail :/")
  }

  # 3. Create directories
  s2ids <- basename(files_points_only_img_folder)
  fpaths <- sprintf("%s/%s/%s/models", output, basename(point), s2ids)
  lapply(
    X = 1:5,
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
  for (index in 1:5) {
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
  for (index in 1:5) {
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
  for (index in 1:5) {
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
  for (index in 1:5) {
    s2xid <- basename(dirname(dirname(fmaskfiles[index])))
    outf <- sprintf("%s/%s/%s/models/fmask.tif", output, basename(point), s2xid)
    Rreshape(fmaskfiles[index], outf, raster_ref)
  }

  # 9. Create npy
  npoutfs <- rep(NA, 5)
  for (index in 1:5) {
    s2folder <- dirname(fpaths[index])
    s2npyid <- basename(s2folder)
    inf <-sprintf("%s/%s/input", point, s2npyid)
    outf <- sprintf("%s/input.npy", s2folder)
    create_npy_file(inf, outf)
    npoutfs[index] <- outf
  }


  # 10. Thumbnail
  for (index in 1:5) {
    s2folder <- dirname(fpaths[index])
    s2npyid <- basename(s2folder)
    thumbnail_in <-sprintf("%s/%s/thumbnails/thumbnail.tif", point, s2npyid)
    thumbnail_out <- sprintf("%s/thumbnail.tif", s2folder)
    thumbnails_creator(thumbnail_in, thumbnail_out, raster_ref)
  }

  # 11. Create UNET
  for (index in 1:5) {
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



points <- list.files("/media/csaybar/Elements SE/cloudSEN12/nolabel/", full.names = TRUE)
output <- "/media/csaybar/Elements/cloudSEN12/nolabel/"
counter <- 1
for (point in points[1:6000]) {
  new_db_migration_local(point, output)
  print(counter)
  counter <- counter + 1
}

