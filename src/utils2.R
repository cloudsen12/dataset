# Functions
cloudsen12_up_local <- function(point_dir, fmask_folder) {
  message("Processing: ", basename(point_dir))

  # 1. Identify the point :)
  point_dir_full_files <- list.files(point_dir, full.names = TRUE)

  # 2. get only the folder with s2/s1 data :D
  files_points_only_img_folder <- point_dir_full_files[grepl("^[0-9]", basename(point_dir_full_files))]
  if (length(files_points_only_img_folder) != 5) {
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
  for (index in 1:5) {
    fullname_s2f <- files_points_only_img_folder[index]
    q60 <- QA60_creator(basename(fullname_s2f), refraster)
    q60_out <- sprintf("%s/target/qa60.tif", fullname_s2f)
    writeRaster(x = q60, filename = q60_out, overwrite = TRUE)
  }

  # 5. Fmask generation
  for (index in 1:5) {
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
