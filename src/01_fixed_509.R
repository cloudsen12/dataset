# Script to fix problems in cloudsen12
library(reticulate)
library(foreach)
library(raster)
library(dplyr)
library(stars)
library(rgee)
library(doMC)

setwd("/home/csaybar/Desktop/extra_lab/unet/")
source_python("main.py")
np <- import("numpy")

registerDoMC(10)
ee_Initialize("gabriela")

local_cloudsen2_metadata <- read.csv("/home/csaybar/Desktop/extra_lab/q60/metadata_s2.csv")
dataset <- "/media/csaybar/Elements SE/cloudSEN12_f/high/"

#1. Check if sen2cor and s2cloudness are the same
check_01  <- function(point) {
  files <- list.files(sprintf("%s/%s", dataset, point), recursive = TRUE, full.names = TRUE)
  r1 <- read_stars(files[grepl("s2cloudless", files)])
  r2 <- read_stars(files[grepl("sen2cor", files)])
  all(all(dim(r1) == dim(r2)), all(st_crs(r1) == st_crs(r2)))
}

#2. Check thumbnail dimension
check_02  <- function(point) {
  files <- list.files(sprintf("%s/%s", dataset, point), recursive = TRUE, full.names = TRUE)
  r1 <- lapply(files[grepl("thumbnail", files)], read_stars)
  lapply(r1, dim)
}

#3. Check fmask4 is consistent with sen2cor
check_03  <- function(point) {
  files <- list.files(sprintf("%s/%s", dataset, point), recursive = TRUE, full.names = TRUE)
  r1 <- read_stars(files[grepl("fmask", files)])
  r2 <- read_stars(files[grepl("sen2cor", files)])
  all(all(dim(r1) == dim(r2)), all(st_crs(r1) == st_crs(r2)))
}

#3. Check qa60 is consistent with sen2cor
check_04  <- function(point) {
  files <- list.files(sprintf("%s/%s", dataset, point), recursive = TRUE, full.names = TRUE)
  r1 <- read_stars(files[grepl("qa60", files)])
  r2 <- read_stars(files[grepl("sen2cor", files)])
  all(all(dim(r1) == dim(r2)), all(st_crs(r1) == st_crs(r2)))
}

#4. Upgrade q60 -----------------------------------
QA60_creator <- function(local_cloudsen2_metadata, point) {
  list_rasters <- list()
  local_cloudsen2_metadata_l <- local_cloudsen2_metadata[local_cloudsen2_metadata$name %in% point,]
  fdir <- sprintf("%s/%s", dataset, local_cloudsen2_metadata_l$name)
  allfiles <- list.files(fdir, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  ref <- allfiles[grepl("sen2cor", allfiles)]
  rr_Ref <- raster(ref[1])
  cente_sp <- st_sfc(st_point(apply(coordinates(rr_Ref), 2, mean)), crs = crs(rr_Ref)@projargs)
  crs_sp <- st_crs(read_stars(ref[1], proxy = TRUE))$epsg

  for (index2 in 1:5) {
    # Sentinel-2 Level 1C
    s2_img <- ee$Image(sprintf("COPERNICUS/S2/%s", local_cloudsen2_metadata_l[, index2 + 1]))$select("QA60")

    # 3. Create a st_point representing the center of the tile (255x255)
    #crs_kernel <- ee$Image(s2_img)$select(0)$projection()$getInfo()$crs
    crs_kernel <- sprintf("EPSG:%s", crs_sp)
    #st_point <- st_sfc(geometry = st_point(c(local_cloudsen2_metadata_l$x, local_cloudsen2_metadata_l$y)), crs = 4326)
    #point_utm <- st_transform(st_point, crs_kernel)
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

    list_rasters[[index2]] <- list(
      qa60 = final_stack,
      point = local_cloudsen2_metadata_l[, 1],
      s2id = local_cloudsen2_metadata_l[, index2 + 1]
    )
  }
  list_rasters
}

export_results <- function(results_qa60, allp_q60) {
  for (index in 1:5) {
    writeRaster(
      x = results_qa60[[index]]$qa60,
      filename = sprintf("%s/%s/models/qa60.tif", allp_q60, results_qa60[[index]]$s2id),
      overwrite = TRUE
    )
  }
}
# --------------------------------------------------


#5. Fmask upgrade -----------------------------------
fmask_up <- function(point) {
  fmask_folder <- "/home/csaybar/Downloads/FMASK4/"
  allfmaskfiles <-  list.files(fmask_folder, "\\.tif$", recursive = TRUE, full.names = TRUE)
  allfmaskfiles_b <- basename(allfmaskfiles)
  fdir <- sprintf("%s/%s", dataset, point)
  fimg <- list.files(fdir)
  fimg_gr_id <- sapply(1:5,
                       function(x) {
                         id <- ee$Image(sprintf("COPERNICUS/S2/%s", fimg[x]))$get("GRANULE_ID")$getInfo()
                         allfmaskfiles[allfmaskfiles_b %in% sprintf("%s_Fmask4.tif", id)][1]
                       }
  )

  if (length(fimg_gr_id) != 5) {
    stop("no fmask4 img")
  }

  to_del <- sprintf("%s/%s/models/fmask.tif", fdir, fimg)
  # get ref
  # crop to 1001
  allfiles <- list.files(fdir, pattern = "\\.tif", full.names = TRUE, recursive = TRUE)
  ref <- allfiles[grepl("sen2cor", allfiles)]
  ddd <- read_stars(ref[1], proxy = TRUE) # check crs
  crs <- st_crs(ddd)$epsg

  r <- raster(ref[1])
  rextent <- extent(r)
  xmin <- rextent[1]
  xmax <- rextent[2]
  ymin <- rextent[3]
  ymax <- rextent[4]

  tempfiletifs <- sprintf("%s.tif", paste0("/home/csaybar/z_", sprintf("%02d",1:5)))
  foreach (index2 = 1:5) %dopar% {
    system(
      sprintf(
        "gdalwarp '%s' '%s' -overwrite -te %s %s %s %s -tr 10 10 -q -t_srs %s",
        fimg_gr_id[index2], tempfiletifs[index2],
        xmin, ymin,
        xmax, ymax, sprintf("EPSG:%s", crs)
      )
    )
    file.copy(
      from = tempfiletifs[index2],
      to = to_del[index2],
      overwrite = TRUE
    )
  }
}

#6. CROP 1001 THUMBNAIL -----------------------------------
thumbnails_creator <- function(point, type = "scribble") {
  # folder1
  # pointf <- sprintf("%s/%s/%s", "/media/csaybar/58059B472A3AA231/", type, point)
  # ips <- list.files(pointf, full.names = TRUE)
  # to_copy <- sprintf("%s/thumbnails/thumbnail.tif", ips[dir.exists(ips)]) %>% sort()
  #
  # pointf2 <- sprintf("%s/high/%s", "/media/csaybar/Elements SE/cloudSEN12_f", point)
  # ips2 <- list.files(pointf2, full.names = TRUE)
  # to_del <-sprintf("%s/thumbnail.tif", ips2) %>% sort()
  # file.copy(to_copy, to_del, overwrite = TRUE)

  pointf2 <- sprintf("%s/%s", dataset, point)
  ips2 <- list.files(pointf2, full.names = TRUE)
  to_del <-sprintf("%s/thumbnail.tif", ips2) %>% sort()

  # crop to 1001
  allfiles <- list.files(pointf2, pattern = "\\.tif", full.names = TRUE, recursive = TRUE)
  ref <- allfiles[grepl("sen2cor", allfiles)]
  ddd <- read_stars(ref, proxy = TRUE) # check crs
  crs <- st_crs(ddd)$epsg

  r <- raster(ref[1])
  centroid <- st_sfc(st_point(c(mean(xFromCol(r)), mean(yFromRow(r)))), crs = crs(r)@projargs)
  ref_b <- st_buffer(centroid, 500*10, endCapStyle="SQUARE")
  xmin <- as.numeric(st_bbox(ref_b))[1] - 5
  ymin <- as.numeric(st_bbox(ref_b))[2] - 5
  xmax <- as.numeric(st_bbox(ref_b))[3] + 5
  ymax <- as.numeric(st_bbox(ref_b))[4] + 5

  tempfiletifs <- sprintf("%s.tif", paste0("/home/csaybar/t_", sprintf("%02d",1:5)))
  foreach (index2 = 1:5) %dopar% {
    system(
      sprintf(
        "gdalwarp '%s' '%s' -overwrite -te %s %s %s %s -tr 10 10 -q -t_srs %s",
        to_del[index2], tempfiletifs[index2],
        xmin, ymin,
        xmax, ymax, sprintf("EPSG:%s", crs)
      )
    )
    file.copy(
      from = tempfiletifs[index2],
      to = to_del[index2],
      overwrite = TRUE
    )
  }
  read_stars(to_del, proxy = TRUE)
}

#7. GLOBAL MODEL CHECK -----------------------------------
check_05 <- function(point) {
  files <- list.files(sprintf("%s/%s", dataset, point), recursive = TRUE, full.names = TRUE)
  r1 <- read_stars(files[grepl("thumbnail", files)])
  r2 <- read_stars(files[!grepl("thumbnail|UNET_|input\\.npy", files)])
  TRUE
}

#8. FIX QA -----------------------------------
fix_01 <- function(point) {
  files <- list.files(sprintf("%s/%s", dataset, point), recursive = TRUE, full.names = TRUE)
  r1 <- stack(files[grepl("qa60", files)])
  r2 <- stack(files[grepl("sen2cor", files)])
  r3 <- resample(r1, r2)
  for (index2 in 1:5) {
    writeRaster(r3[[index2]], files[grepl("qa60", files)][index2], overwrite = TRUE)
  }
}

#9. clic 509 -----------------------------------
eliminate_one_pixel <- function(point) {
  files <- list.files(sprintf("%s/%s", dataset, point), recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("thumbnail|input\\.npy|UNET_", files)]
  stk <- stack(files)
  xmin <- extent(stk)[1] + 10
  xmax <- extent(stk)[2] - 10
  ymin <- extent(stk)[3] + 10
  ymax <- extent(stk)[4] - 10
  rr <- raster(
    nrows = 509, ncols = 509,
    xmn = xmin, xmx = xmax,
    ymn = ymin, ymx = ymax,
    crs = crs(r)
  )
  for (index in 1:20) {
    r <- stk[[index]]
    rr[] <- r[2:510,2:510]
    writeRaster(rr, files[index], overwrite = TRUE)
  }
}

# 10. clic numpy 509 -----------------------------------
eliminate_one_pixel_np <- function(file) {
  files <- list.files(sprintf("%s/%s", dataset, file), recursive = TRUE, full.names = TRUE)
  npys <- files[grepl("input\\.npy$", files)]
  for (index2 in 1:5) {
    s2cube <- np$load(npys[index2])[, 2:510, 2:510]
    np$save(npys[index2], s2cube)
  }
}

# 11. check np  ---------------------------------
check_06 <- function(point) {
  files <- list.files(sprintf("%s/%s", dataset, point), recursive = TRUE, full.names = TRUE)
  npys <- files[grepl("input\\.npy$", files)]
  list_df <- list()
  for (index2 in 1:5) {
    list_df[[index2]] <- dim(np$load(npys[index2]))
  }
  do.call(rbind, list_df) %>% as.data.frame()
}


# 12. From np to stack -----------------------------------
from_npy_to_stack <- function(point) {
  files <- list.files(sprintf("%s/%s", dataset, point), recursive = TRUE, full.names = TRUE)
  npys <- files[grepl("input\\.npy$", files)]
  rbase <- raster(files[grepl("sen2cor", files)][1])

  rrlist <- list()
  for (index2 in 1:5) {
    nparray <- np$load(npys[index2])
    rlist <- list()
    for (index3 in 1:20) {
      rbase[] <- t(nparray[index3,,])
      rlist[[index3]] <- rbase
    }
    rlistf <- stack(rlist)
    names(rlistf) <- c(
      "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11",
      "B12", "VH", "VV", "angle", "CDI", "cloudshadow_direction", "elevation",
      "landuse"
    )
    rrlist[[index2]] <- rlistf
  }
  rrlist
}

# 13. unet  -----------------------------------
unet_gen <- function(point) {
  files <- list.files(sprintf("%s/%s", dataset, point), recursive = TRUE, full.names = TRUE)
  npys <- files[grepl("input\\.npy$", files)]
  sen2cor_input <- files[grepl("sen2cor", files)][1]
  base <- raster(sen2cor_input)
  for (index2 in 1:5) {
    npy_file_input <- npys[index2]
    model_results01 <- py$unet_main(npy_file = npy_file_input, namemodel = "rgbiswir")
    model_results02 <- py$unet_main(npy_file = npy_file_input, namemodel = "rgbi")

    model_results_raster_01 <- base
    model_results_raster_01[] <- model_results01

    model_results_raster_02 <- base
    model_results_raster_02[] <- model_results02

    writeRaster(
      x = model_results_raster_01,
      filename = sprintf("%s/models/UNET_rgbiswir.tif", dirname(npy_file_input)),
      overwrite = TRUE
    )
    writeRaster(
      x = model_results_raster_02,
      filename = sprintf("%s/models/UNET_rgbi.tif", dirname(npy_file_input)),
      overwrite = TRUE
    )
  }
}

#########################
#RUN
#########################
all_points <- list.files(dataset)
results <- foreach(index = 1:192) %do% {
  print(index)
#ddd <- list()
#for (index in 1:192) {
  # check_01(point = all_points[index]) # # Check if sen2cor and s2cloudness same area.
  # check_02(all_points[index]) # Check thumbnail dimension
  # check_03(point = all_points[index]) # Check fmask4 is consistent with sen2cor
  # check_04(point = all_points[index]) # Check qa60 is consistent with sen2cor

  # QA60 ----------------------------------------
  #print(index)
  #results_qa60 <- QA60_creator(local_cloudsen2_metadata, all_points[index], "scribble")
  #allp_q60 <- sprintf("/media/csaybar/58059B472A3AA231/scribble/%s", results_qa60[[1]]$point)
  #export_results(results_qa60, allp_q60)
  # ---------------------------------------------

  # FMASK ----------------------------------------
  # fmask_up(point = all_points[index], type = "scribble")

  # THUMBNAIL -----------------------------------
  # print(index)
  #thumbnails_creator(all_points[index], type = "scribble")

  # GLOBAL CHECK -----------------------------------
  # check_05(point = all_points[index])
  # fix_01(all_points[index]) # FIX QA

  # CLIC  509X509-----------------------------------
  eliminate_one_pixel(point = all_points[index])

  # CLIC np  509X509 -------------------------------
  # print(index)
  # ddd[[index]] <- check_06(point = all_points[index])
  # eliminate_one_pixel_np(file = all_points[index])

  # RUN UNET  -------------------------------
  # print(index)
  # unet_gen(point = all_points[index])
}

sum(unlist(results))
# final_r <- do.call(rbind, lapply(1:1808, function(x) do.call(rbind, results[[x]]))) %>% as.data.frame()
# table(final_r$y)
# final_r <- do.call(rbind, ddd)

# all_points <- list.files(dataset, full.names = TRUE)
# all_bad <- list()
# for (index in 1:1808) {
#   print(index)
#   files <- list.files(all_points[index], pattern = "\\.tif.aux.xml$", recursive = TRUE, full.names = TRUE)
#   all_bad[[index]] <- files
# }
#unlist(all_bad)
