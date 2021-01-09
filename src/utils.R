#' CloudSEN12 thumbnail creator
#'
#' Create several thumbnail to find the clouds closest to the theoretical probabilities
#' All Sentinel2 images have a Sentinel1 pair with no more than 2.5 days of delay.
#'
#' @param n_images  Numeric. Number of images to download, if n_images = "max" the
#' function will download all the images available.
#' @param kernel_size Size of the kernel.
#' @param data_range Range of dates to obtain images.
#' @param output Folder where to save the results.
#'
select_dataset_thumbnail_creator <- function(cloudsen2_row,
                                             n_images = "max",
                                             kernel_size = c(255, 255),
                                             data_range = c("2018-01-01", "2020-07-31"),
                                             output = "results/") {
  # 1. Create output directory
  dir.create(output, showWarnings = FALSE)

  # 2. Create a point which represent the center of the chip (from local to Earth Engine)
  point <- ee$Geometry$Point(cloudsen2_row$geometry[[1]])

  # 3. Create a S2 ImageCollection.
  s2Sr <- ee$ImageCollection("COPERNICUS/S2_SR") %>%
    ee$ImageCollection$filterBounds(point) %>%
    ee$ImageCollection$filterDate(data_range[1], data_range[2])

  # 4. Create a S1 ImageCollection.
  data_range2 <- c(as.Date(data_range[1]) - 5, as.Date(data_range[2]) + 5) %>% as.character()
  s1_grd <- ee$ImageCollection("COPERNICUS/S1_GRD") %>%
    ee$ImageCollection$filterBounds(point) %>%
    ee$ImageCollection$filterDate(data_range2[1], data_range2[2]) %>%
    # Filter to get images with VV and VH dual polarization.
    ee$ImageCollection$filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV")) %>%
    ee$ImageCollection$filter(ee$Filter$listContains('transmitterReceiverPolarisation', "VH")) %>%
    # Filter to get images collected in interferometric wide swath mode.
    ee$ImageCollection$filter(ee$Filter$eq("instrumentMode", "IW"))

  # 5. Get dates from all the images which are part of S1 and S2 ImageCollection.
  s2_dates <- ee_get_date_ic(s2Sr) %>% as_tibble()
  s1_dates <- ee_get_date_ic(s1_grd) %>% as_tibble()
  sx_fx <- function(x) min(abs(s2_dates$time_start[x] - s1_dates$time_start))
  mindays <- sapply(seq_len(nrow(s2_dates)), sx_fx)
  valid_s2 <- s2_dates[mindays < 2.5,] # Pick up images with no more than 2.5 days of delay

  # Display the number of available images with an S2/S1 pair
  message(sprintf("Number of images: %s", nrow(valid_s2)))

  if (n_images == "max") {
    n_images <- nrow(valid_s2)
  }

  # 6. Get the CRS of this specific point (5 images)
  img_crs <- s2Sr$first()$select(0)$projection()$getInfo()[["crs"]]

  # 7. shuffle valid images
  images_position <- sample(nrow(valid_s2), nrow(valid_s2))

  # if (length(images_position) < 50) {
  #   warning("Insufficient number of images ... PLEASE REPORT!")
  # }

  # 8. Create a folder to save results
  dir_id <- sprintf("%s/point_%04d",output, cloudsen2_row$id)
  dir.create(dir_id, showWarnings = FALSE)

  # 9. Download all the image thumbnails
  counter <- 0
  for (r_index in images_position) {
    if (counter == n_images) {
      break
    }
    # 9.1 Select the image
    img_to_download <- ee$Image(valid_s2$id[r_index])
    img_id <- basename(valid_s2$id[r_index])

    # 9.2 Download the S2 thumbnail image (as a json)
    s2_img_array <- img_to_download %>%
      ee$Image$select(c("B4", "B3", "B2")) %>%
      ee$Image$addBands(ee$Image$pixelCoordinates(projection = img_crs)) %>%
      # ee$Image$reproject(crs = img_crs, scale = 10) %>%
      ee$Image$neighborhoodToArray(
        kernel = ee$Kernel$rectangle(kernel_size[1], kernel_size[2], "pixels")
      ) %>%
      ee$Image$sampleRegions(ee$FeatureCollection(point),
                             projection = img_crs,
                             scale = 10) %>%
      ee$FeatureCollection$getInfo()


    # Some image have FULL NA values if that occurs skip
    band_names <- try(
      expr = names(s2_img_array$features[[1]]$properties),
      silent = TRUE
    )
    if (class(band_names) == "try-error") {
      next
    }

    message(
      sprintf("Processing point [%s] image [%s] ... please wait",
              cloudsen2_row$id, counter)
    )

    # 9.3 From list to data_frame
    extract_fn <- function(x) as.numeric(unlist(s2_img_array$features[[1]]$properties[x]))
    image_as_df <- do.call(cbind,lapply(band_names, extract_fn))
    colnames(image_as_df) <- band_names
    image_as_tibble <- as_tibble(image_as_df)

    # If all values in the image are zero, skip
    if (sum(image_as_tibble[["B2"]] == 0) > 0) {
      next
    } else {
      counter <- counter + 1
    }

    # 9.4 From data_frame to sp; From sp to raster
    coordinates(image_as_tibble) <- ~x+y
    sf_to_stack <- function(x) rasterFromXYZ(image_as_tibble[x])
    final_stack <- stack(lapply(names(image_as_tibble), sf_to_stack))

    # 9.4 Create a plot and download
    crs(final_stack) <- st_crs(img_crs)$proj4string
    png(sprintf("%s/%s.png", dir_id, img_id), 1000, 1000)
    max_value <- max(maxValue(final_stack))
    plotRGB(final_stack/max_value, r = 3, g = 2, b = 1, scale = 1)
    dev.off()
  }

  # 10. Create the metadata.json file
  metadata_dataset_creator(
    cloudsen2_row = cloudsen2_row,
    output = output
  )
}


#' Function to create metadata  (i.e. metadata_1500.json)
#'
#' This function is used in select_dataset_thumbnail_creator to create the
#' final json file to be fill out for the labelers.
#'
#' @param jsonfile metadata (*.json) with the sentinel2 ID.
#' @param kernel_size Size of the kernel.
#' @param output_final Folder where to save the results.
#'
dataset_creator_chips <- function(jsonfile,
                                  kernel_size = c(255, 255),
                                  output_final = "cloudsen12/") {

  point_name <- paste0("point_", gsub("[a-zA-Z]|_|\\.","", basename(jsonfile)))
  # 1. Read JSON file
  jsonfile_r <- jsonlite::read_json(jsonfile)

  # 2. Identify all the S2 images
  s2_ids <- sprintf("COPERNICUS/S2/%s", names(jsonfile_r)[1:5])

  # 3. Create a st_point which represent the center of the chip
  st_point <- st_sfc(geometry = st_point(c(jsonfile_r$x, jsonfile_r$y)), crs = 4326)
  crs_kernel <- ee$Image(s2_ids[1])$select(0)$projection()$getInfo()$crs
  point_utm <- st_transform(st_point, crs_kernel)
  ee_point <- ee$Geometry$Point(point_utm[[1]], proj = crs_kernel)

  # 4. Download each image for the specified point
  for (s2_id in s2_ids) {
    message(sprintf("Downloading: %s", s2_id))

    # 3.1 S2 ID and dates
    s2_img <- ee$Image(s2_id)
    s2_date <- ee_get_date_img(s2_img)[["time_start"]]

    # 3.2 S1 ID
    s1_id <- ee_get_s1(point = ee_point, s2_date = s2_date)
    s1_img <- ee$Image(s1_id)

    # 3.3 Create an Image collection with S2, S1 and cloud mask information
    s2_fullinfo <- ee_merge_s2_full(s2_id, s1_id, s2_date)


    # 3.4 Create a 511x511 tile (list -> data_frame -> sp -> raster)
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
    extract_fn <- function(x) as.numeric(unlist(s2_img_array$features[[1]]$properties[x]))
    image_as_df <- do.call(cbind,lapply(band_names, extract_fn))
    colnames(image_as_df) <- band_names
    image_as_tibble <- as_tibble(image_as_df)
    coordinates(image_as_tibble) <- ~x+y
    sf_to_stack <- function(x) rasterFromXYZ(image_as_tibble[x])
    final_stack <- stack(lapply(names(image_as_tibble), sf_to_stack))
    crs(final_stack) <- st_crs(crs_kernel)$proj4string

    # 3.5 Prepare folders for iris
    output_final_d <- sprintf("%s/dataset", output_final)
    output_final_folder <- sprintf("%s/dataset/%s/%s", output_final, point_name, basename(s2_id))

    metadata_main <- sprintf("%s/cloud_segmentation_%s.json", output_final, point_name)
    metadata_spec <- sprintf("%s/dataset/%s/%s/metadata.json", output_final, point_name, basename(s2_id))

    dir.create(sprintf("%s/input", output_final_folder), showWarnings = FALSE, recursive = TRUE)
    dir.create(sprintf("%s/target", output_final_folder), showWarnings = FALSE, recursive = TRUE)
    dir.create(sprintf("%s/thumbnails", output_final_folder), showWarnings = FALSE, recursive = TRUE)

    # 3.6 Save all features inside the input folder
    bandnames <- c(paste0("B",1:8), "B8A", paste0("B", 9:12), "CDI", "VV", "VH", "angle", "elevation", "landuse")
    input_data <- raster::stack(
      final_stack[[1:13]]/10000, final_stack[[14]], final_stack[[15:17]], final_stack[[22:23]]
    )
    input_spec <- sprintf("%s/input/%s.tif", output_final_folder, bandnames)
    lapply(1:19, function(x) writeRaster(input_data[[x]], input_spec[x], overwrite = TRUE))

    # 3.7 Save all cloud mask inside the target folder
    # 18-19 -> cmask_s2cloudness| cmask_s2cloudness_reclass (0,1)
    # 20-21 -> cmask_sen2cor | cmask_sen2cor_reclass (0,1,2)
    bandnames <- c("s2cloudness", "s2cloudness_reclass", "sen2cor", "sen2cor_reclass")
    benchmarch_data <- final_stack[[18:21]]
    target_spec <- sprintf("%s/target/%s.tif", output_final_folder, bandnames)
    lapply(1:4, function(x) writeRaster(benchmarch_data[[x]], target_spec[x], overwrite = TRUE))

    # 3.8 Create cloud-segmentation.json (main file in iris software)
    ee_create_cloudseg(path = metadata_main)

    # 3.9 Create metadata.json for each file
    ee_create_metadata(
      id = basename(s2_id),
      point = c(jsonfile_r$y, jsonfile_r$x),
      path = metadata_spec
    )
  }

  # 4 Save geometry
  roi <- extent(final_stack[[1]]) %>%
    st_bbox() %>%
    st_as_sfc()
  st_crs(roi) <- crs_kernel
  write_sf(roi, sprintf("%s/%s.gpkg", dirname(output_final_folder), point_name))
}

#' Function to create metadata  (i.e. metadata_1500.json)
#'
#' This function is used in select_dataset_thumbnail_creator to create the
#' final json file to be fill out for the labelers.
#'
#' @noRd
metadata_dataset_creator <- function(cloudsen2_row, output) {
  dir_name_point <- sprintf("%s/point_%04d/", output, cloudsen2_row$id)
  cloudsen2_row_no_sf <- st_drop_geometry(cloudsen2_row)
  cprob_n <- names(cloudsen2_row_no_sf[4:8])
  pprob_v <- as.numeric(cloudsen2_row_no_sf[4:8])

  row_id <- cloudsen2_row$id
  metadata_point <- local_cloudsen2_points[row_id,]
  coord_xy <- as.numeric(metadata_point$geometry[[1]])
  param_id <- list()
  for (index in seq_along(cprob_n)) {
    param_id[[index]] <- list(
      cloud_type = NA,
      cloud_height = NA,
      cloud_thickness = NA,
      potential_cloud_coverage = pprob_v[index]
    )
  }
  names(param_id) <- cprob_n
  param_id$surface_type <- as.character(metadata_point$type)
  param_id$x <- coord_xy[1]
  param_id$y <- coord_xy[2]
  param_id$comments <- "PUT_HERE_YOUR_COMMENT"
  jsonlite::write_json(
    x = param_id,
    path = sprintf("%s/metadata_%04d.json", dir_name_point, cloudsen2_row$id),
    pretty = TRUE,
    auto_unbox = TRUE
  )
}

#' Get the SENTINEL1 image for a specific SENTINEL2
#'
#' This function is used in select_dataset_thumbnail_creator to create the
#' final json file to be fill out for the labelers.
#'
#' @noRd
ee_get_s1 <- function(point, s2_date, range = 2.5) {

  # 1. Defining temporal filter
  s1_date_search <- list(
    init_date = (s2_date - lubridate::hours(range * 24)) %>% rdate_to_eedate(),
    last_date = (s2_date + lubridate::hours(range * 24)) %>% rdate_to_eedate()
  )

  # 2. Load S1 data
  s1_grd <- ee$ImageCollection("COPERNICUS/S1_GRD") %>%
    ee$ImageCollection$filterBounds(point) %>%
    ee$ImageCollection$filterDate(s1_date_search[[1]], s1_date_search[[2]]) %>%
    # Filter to get images with VV and VH dual polarization.
    ee$ImageCollection$filter(ee$Filter$listContains("transmitterReceiverPolarisation", "VV")) %>%
    ee$ImageCollection$filter(ee$Filter$listContains('transmitterReceiverPolarisation', "VH")) %>%
    # Filter to get images collected in interferometric wide swath mode.
    ee$ImageCollection$filter(ee$Filter$eq("instrumentMode", "IW"))

  # 3. get dates and ID
  s1_grd_id <- tryCatch(
    expr = ee_get_date_ic(s1_grd),
    error = function(e) {
      message(
        "Not found any Sentinel-1 within 2.5 days ...,",
        "This is so weird :| ... Report to Cesar :)"
      )
      stop(e)
    }
  )

  # 4. Get the nearest image
  row_position <- which.min(abs(s1_grd_id$time_start - s2_date))
  s1_grd_id[row_position,][["id"]]
}


#' Merge SENTINEL2, SENTINEL1 and cloud label (sen2cor, s2cloudness) predictions
#'
#' @noRd
ee_merge_s2_full <- function(s2_id, s1_id, s2_date) {

  # 1. Define all the images that will be inside input.tif (ImageCollection)
  ## S1
  s1_grd <- ee$Image(s1_id)

  ## S2-level 2A (sen2cor)
  s2_2a <- ee$Image(sprintf("COPERNICUS/S2_SR/%s", basename(s2_id)))

  ## S2-level 1C
  s2_1c <- ee$Image(sprintf("COPERNICUS/S2/%s", basename(s2_id)))
  ### Estimate CDI
  s2_cdi <- ee$Algorithms$Sentinel2$CDI(s2_1c)

  ## DEM data (MERIT)
  extra_dem <- cloudsen12_dem()

  ## LandUSE data (copernicus)
  extra_LC <- cloudsen12_lc()

  # 2. Define all the images that will be inside target.tif (ImageCollection)
  ## Cloud mask according to Zupanc et al. 2019
  s2_cloud <- ee$Image(sprintf("COPERNICUS/S2_CLOUD_PROBABILITY/%s", basename(s2_id)))
  boxcar1 <- ee$Kernel$square(radius = 4, units = 'pixels')
  boxcar2 <- ee$Kernel$square(radius = 2, units = 'pixels')
  s2cloudness_prob <- s2_cloud$rename("cmask_s2cloudness")
  s2cloudness_prob_reclass <- s2_cloud %>%
    ee$Image$convolve(boxcar1) %>%
    ee$Image$focal_max(kernel = boxcar2, iterations = 1) %>%
    ee$Image$gte(70) %>%
    ee$Image$rename("cmask_s2cloudness_reclass")

  # 3. Reclass sen2cor (SLC, sentinel2 level2A)
  ## 4,5,6,11 -> clear
  ## 8,9,10 -> cloud
  ## 2, 3 -> cloud shadows
  ## 1, 7 -> no data
  ## OBS: In Earth Engine "no data" values are masked out.
  s2_scl <- s2_2a$select("SCL")$rename("cmask_sen2cor")
  s2_scl_reclass <- s2_2a$select("SCL")$remap(
    c(4, 5, 6, 11, 8, 9, 10, 2, 3, 1, 7),
    c(0, 0, 0, 0, 1, 1, 1, 2, 2, 1, 1)
  )$rename("cmask_sen2cor_reclass")

  # 5. Merge S2_1C (B.*) + CDI +  S2_2A (SLC) and S2_CLOUD_PROBABILITY (cloud_prob)
  s2_1c %>%
    ee$Image$select("B.*") %>%
    ee$Image$addBands(s2_cdi) %>%
    ee$Image$addBands(s1_grd) %>%
    ee$Image$addBands(s2cloudness_prob) %>%
    ee$Image$addBands(s2cloudness_prob_reclass) %>%
    ee$Image$addBands(s2_scl) %>%
    ee$Image$addBands(s2_scl_reclass) %>%
    ee$Image$addBands(extra_dem) %>%
    ee$Image$addBands(extra_LC)
}

#' Create IRIS main json
#' @noRd
#'
ee_create_cloudseg <- function(path) {
  point_name <- gsub("\\.json$","",paste0(strsplit(basename(path), "_")[[1]][3:4],collapse = "_"))
  cseg_list <- list(
    name = point_name,
    authentication_required = TRUE,
    images = list(
      path = list(
        B1 = sprintf("dataset/%s/{id}/input/B1.tif", point_name),
        B2 = sprintf("dataset/%s/{id}/input/B2.tif", point_name),
        B3 = sprintf("dataset/%s/{id}/input/B3.tif", point_name),
        B4 = sprintf("dataset/%s/{id}/input/B4.tif", point_name),
        B5 = sprintf("dataset/%s/{id}/input/B5.tif", point_name),
        B6 = sprintf("dataset/%s/{id}/input/B6.tif", point_name),
        B7 = sprintf("dataset/%s/{id}/input/B7.tif", point_name),
        B8 = sprintf("dataset/%s/{id}/input/B8.tif", point_name),
        B8A = sprintf("dataset/%s/{id}/input/B8A.tif", point_name),
        B9 = sprintf("dataset/%s/{id}/input/B9.tif", point_name),
        B10 = sprintf("dataset/%s/{id}/input/B10.tif", point_name),
        B11 = sprintf("dataset/%s/{id}/input/B11.tif", point_name),
        B12 = sprintf("dataset/%s/{id}/input/B12.tif", point_name),
        CDI = sprintf("dataset/%s/{id}/input/CDI.tif", point_name),
        VV = sprintf("dataset/%s/{id}/input/VV.tif", point_name),
        VH = sprintf("dataset/%s/{id}/input/VH.tif", point_name),
        angle = sprintf("dataset/%s/{id}/input/angle.tif", point_name),
        elevation = sprintf("dataset/%s/{id}/input/elevation.tif", point_name),
        landuse = sprintf("dataset/%s/{id}/input/landuse.tif", point_name),
        s2cloudness = sprintf("dataset/%s/{id}/target/s2cloudness.tif", point_name),
        s2cloudness_reclass = sprintf("dataset/%s/{id}/target/s2cloudness_reclass.tif", point_name),
        sen2cor = sprintf("dataset/%s/{id}/target/sen2cor.tif", point_name),
        sen2cor_reclass = sprintf("dataset/%s/{id}/target/sen2cor_reclass.tif", point_name)
      ),
      shape = c(511,511),
      metadata = sprintf("dataset/%s/{id}/metadata.json", point_name)
    ),
    segmentation = list(
      path = sprintf("dataset/%s/{id}/target/manual.tif", point_name),
      mask_encoding = "rgb",
      mask_area = c(0, 0, 511, 511),
      score = "f1",
      pending_threshold = 1,
      test_images = NA
    ),
    classes = list(
      list(
        name = "Clear",
        description = "All clear pixels, i.e. without cloud contamination or cloud shadows.",
        colour = c(255,255,255,0),
        user_colour = c(0,255,255,70)
      ),
      list(
        name = "Thick Cloud",
        description = "All cloudy pixels covered by thick clouds (does not include semi-transparent clouds or cloud shadows).",
        colour = c(255, 255, 0, 70)
      ),
      list(
        name = "Thin Cloud",
        description = "Clouds that are semi-transparent, i.e. one can see land or sea surfaces through them. If a thin cloud lays over a thick cloud, please paint them with the <i>Thick Cloud</i> class.",
        colour = c(0, 255, 0, 70)
      ),
      list(
        name = "Cloud Shadows",
        description = "All pixels contaminated by cloud shadows (not terrain shadows).",
        colour = c(255, 0, 0, 70)
      ),
      list(
        name = "No data",
        description = "Reserved for no data pixels, e.g. pixels outside of the satellite's swath.",
        colour = c(50, 50, 255, 70)
      )
    ),
    views = list(
      Cirrus = list(
        description = "Cirrus and high clouds are red.",
        type = "image",
        data = "$B11.B1**0.8*5",
        cmap = "jet"
      ),
      cloud_index = list(
        description = "Cloud Displacement Index, clouds are red.",
        type = "image",
        data = "$CDI.B1*-1"
      ),
      RGB = list(
        description = "Normal RGB image.",
        type = "image",
        data = c("$B5.B1", "$B3.B1", "$B2.B1")
      ),
      NRGB = list(
        description = "Near-Infrared RGB image.",
        type = "image",
        data = c("$B5.B1*1.5", "$B3.B1*1.5", "$B2.B1*1.5")
      ),
      Snow = list(
        description = "Small ice crystals in high-level clouds appear reddish-orange or peach, and thick ice snow looks vivid red (or red-orange). Bare soil appears bright cyan and vegetation seem greenish in the image. Water on the ground is very dark as it absorbs the SWIR and the red, but small (liquid) water drops in the clouds scatter the light equally in both visible and the SWIR, and therefore it appears white. Water Sediments are displayed as dark red.",
        type = "image",
        data = c("$B10.B1", "$B11.B1", "$B12.B1")
      ),
      "Sentinel-1" = list(
        description = "RGB of VH, VV and VH-VV.",
        type = "image",
        data = c("$VH.B1", "$VV.B1", "$VH.B1-$VV.B1")
      ),
      Bing = list(
        description = "Aerial Imagery",
        type = "bingmap"
      ),
      elevation = list(
        description = "Elevation values",
        type = "image",
        data = "$elevation.B1"
      ),
      NDVI = list(
        description = "NDVI values",
        type = "image",
        data = "($B5.B1 - $B4.B1)/($B5.B1 + $B4.B1)"
      ),
      sen2cloudness = list(
        description = "Sen2Cloudness Probability",
        type = "image",
        data = "$s2cloudness.B1"
      ),
      sen2cloudness_reclass = list(
        description = "Sen2Cloudness Probability Reclass (BLUE-->CLEAR; RED -> CLOUD)",
        type = "image",
        data = "$s2cloudness_reclass.B1"
      ),
      sen2cor = list(
        "description" = "sen2cor classes (BLUE--> CLEAR; GREEN -> CLOUD; RED -> CLOUD SHADOW)",
        "type" = "image",
        "data" = "$sen2cor_reclass.B1"
      )
    ),
    view_groups = list(
      default = c("Cirrus", "RGB", "Snow")
    )
  )
  jsonlite::write_json(
    x = cseg_list,
    path = path,
    pretty = TRUE,
    auto_unbox = TRUE
  )
}

#' Create IRIS relative json (for each image)
#' @noRd
#'
ee_create_metadata <- function(id, point, path) {
  scene_list <- list(
    spacecraft_id = "Sentinel2/Sentinel1",
    scene_id = id,
    location = point,
    resolution = 10.0
  )
  jsonlite::write_json(
    x = scene_list,
    path = path,
    pretty = TRUE,
    auto_unbox = TRUE
  )
}


#' Funcion de pato para subir a drive no es mia :X
#' @noRd
#'
drive_upload_full <- function(from, to) {
  drive_mkdir(name = basename(from), path = to)
  map(
    list.files(from, full.names = T, recursive = F),
    ~ drive_upload(
      .x, verbose = FALSE,
      path = sprintf("%1s%2s/", to, basename(from))
    )
  )
}

#' Merge MERIT + CrySat2 DEM
#' @noRd
#'
cloudsen12_dem <- function() {
  merit_dem <- ee$Image("MERIT/Hydro/v1_0_1")$select("elv")
  cryosat2_dem <- ee$Image("CPOM/CryoSat2/ANTARCTICA_DEM")$select("elevation")
  cloudsen12_dem <- ee$Image$add(merit_dem$unmask(0), cryosat2_dem$unmask(0))
  cloudsen12_dem
}


#' Proba-V-C3 + ANTARCTICA (70 -> snow)
#' @noRd
#'
cloudsen12_lc <- function() {
  extra_LC <- ee$ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global") %>%
    ee$ImageCollection$filterDate("2018-12-31", "2019-01-02") %>%
    ee$ImageCollection$first() %>%
    ee$Image$select("discrete_classification") %>%
    ee$Image$rename("land_cover")
  cryosat2_dem <- ee$Image("CPOM/CryoSat2/ANTARCTICA_DEM")$select("elevation")
  antarctica_LC <- cryosat2_dem$multiply(0)$add(70)$unmask(0, sameFootprint = FALSE)
  ee$Image$add(extra_LC$unmask(0), antarctica_LC)
}

# Search metajson in roy folder :)
search_metajson <- function(pattern, clean = TRUE) {
  drive <- sprintf("%s/drive_dataset.Rdata", tempdir())

  if (clean) {
    suppressWarnings(file.remove(drive))
  }

  if (!file.exists(drive)) {
    # 5. List all the metadata
    drive_jsonfile <- drive_ls(
      path = as_id("1fBGAjZkjPEpPr0p7c-LtJmfbLq3s87RK")
    )
    save(drive_jsonfile, file = drive)
  } else {
    load(drive)
  }
  drive_jsonfile_s <- drive_jsonfile[drive_jsonfile$name %in% pattern, ]
  if (nrow(drive_jsonfile_s) == 0) {
    stop("no metadata found")
  }
  jsonfile <- try(
    drive_download(
      file = drive_jsonfile_s,
      path = paste0(tempdir(),"/", drive_jsonfile_s$name),
      overwrite = TRUE
    )
  )
  stop_5 <- 0
  while((class(jsonfile)[1] == "try-error") | stop_5 == 5) {
    jsonfile <- try(
      drive_download(
        file = drive_jsonfile_s,
        path = paste0(tempdir(),"/", drive_jsonfile_s$name),
        overwrite = TRUE
      )
    )
    stop_5 = stop_5 + 1
  }
  paste0(tempdir(),"/", drive_jsonfile_s$name)
}

