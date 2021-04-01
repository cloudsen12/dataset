#' CloudSEN12 thumbnail creator
#'
#' Create several thumbnail to find the clouds closest to the theoretical probabilities
#' All sentinel2 images have a Sentinel1 pair with no more than 2.5 days of delay.
#'
#' @param n_images  Numeric. Number of images to download, if n_images = "max" the
#' function will download all the available images.
#' @param kernel_size Size of the kernel.
#' @param data_range Range of dates to obtain images.
#' @param output Folder where to save the results.
#' @noRd
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
  s2Sr <- ee$ImageCollection("COPERNICUS/S2") %>%
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
  s1_dates <- tryCatch(
    expr = ee_get_date_ic(s1_grd) %>% as_tibble(),
    error = function(e) stop("Sentinel-1 no images ... Please Report!: ", e)
  )
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


#' Main function to create CLOUDSEN12 data
#'
#' This function is used to download all the images.
#'
#' @param jsonfile metadata (*.json) Extract files from Roy Drive.
#' @param kernel_size Size of the kernel.
#' @param output_final Folder where to save the results.
#'
dataset_creator_chips <- function(jsonfile,
                                  sp_db,
                                  kernel_size = c(255, 255),
                                  output_final = "cloudsen12/") {
  point_name <- paste0("point_", gsub("[a-zA-Z]|_|\\.","", basename(jsonfile)))

  # 1. Read JSON file
  jsonfile_r <- jsonlite::read_json(jsonfile)

  # 2. Identify all the S2 images
  s2_idsposition <- which(sapply(strsplit(names(jsonfile_r), "_"), length) == 3)
  s2_ids <- sprintf("COPERNICUS/S2/%s", names(jsonfile_r)[s2_idsposition])

  # 3. Create a st_point representing the center of the tile (255x255)
  st_point <- st_sfc(geometry = st_point(c(jsonfile_r$x, jsonfile_r$y)), crs = 4326)
  crs_kernel <- ee$Image(s2_ids[1])$select(0)$projection()$getInfo()$crs
  point_utm <- st_transform(st_point, crs_kernel)
  ee_point <- ee$Geometry$Point(point_utm[[1]], proj = crs_kernel)

  # 4. Download each image for the specified point
  counter <- 0
  s2_dates <- rep(NA,5)
  s1_dates <- rep(NA,5)
  s1_ids <- rep(NA,5)
  s1_area_per <- rep(NA,5)
  land_use <- rep(NA,5)
  elevation <- rep(NA,5)
  shadow_dir <- rep(NA,5)

  for (s2_id in s2_ids) {
    message(sprintf("Downloading: %s", s2_id))

    # 4.1 S2 ID and dates
    s2_img <- ee$Image(s2_id)
    s2_date <- ee_get_date_img(s2_img)[["time_start"]]
    s2_dates[counter + 1] <- s2_date %>% as.character()

    # 4.2 S1 ID
    s1_id <- ee_get_s1(point = ee_point, s2_date = s2_date)
    s1_img <- ee$Image(s1_id)
    s1_date <- ee_get_date_img(s1_img)[["time_start"]]
    s1_ids[counter + 1] <- s1_id
    s1_dates[counter + 1] <- s1_date %>% as.character()

    # 4.3 Create an Image collection with S2, S1 and cloud mask information
    s2_s1_img <- ee_merge_s2_full(s2_id, s1_id, s2_date)

    # 4.4 Add the shadow direction
    s2_fullinfo <- s2_s1_img %>%
      ee$Image$addBands(shadow_direction(image = s2_img))

    # 4.5 IPL_UV algorithm ... exist enough images?
    IPL_multitemporal_cloud_logical <- ee_upl_cloud_logical(
      sen2id = basename(s2_id),
      roi =  s2_img$geometry()
    )

    # 4.6 If IPL_multitemporal_cloud_logical is TRUE add to the EE dataset
    if (IPL_multitemporal_cloud_logical) {
      IPL_multitemporal_cloud <- ee_upl_cloud(
        sen2id = basename(s2_id),
        roi =  s2_img$geometry()
      ) %>% ee$Image$unmask(-99, sameFootprint = FALSE)
      s2_fullinfo <- s2_fullinfo %>%
        ee$Image$addBands(IPL_multitemporal_cloud)
    }

    # 4.7 Create a 511x511 tile (list -> data_frame -> sp -> raster)
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

    # 4.9 Estimate dataset_values
    ## s1_area_per
    s1_raster_values <- getValues(final_stack[[16]])
    na_pixels <- sum(is.na(s1_raster_values))
    n99_pixels <- sum(s1_raster_values == -99, na.rm = TRUE)
    s1_area <- 100 - (na_pixels + n99_pixels)/ncell(s1_raster_values)*100
    s1_area_per[counter + 1] <- round(s1_area, 2)
    ## landuse
    land_use[counter + 1] <- Modes(getValues(final_stack[[23]]))
    ## elevation
    elevation[counter + 1] <- round(mean(getValues(final_stack[[22]])), 2)
    ## shadow direction
    shadow_dir[counter + 1] <- round(mean(getValues(final_stack[[24]])), 2)


    # 4.8 Prepare folders for iris
    output_final_d <- sprintf("%s/dataset", output_final)
    output_final_folder <- sprintf("%s/dataset/%s/%s", output_final, point_name, basename(s2_id))

    metadata_main <- sprintf("%s/dataset/%s/cloud_segmentation_%s.json", output_final, point_name, point_name)
    metadata_spec <- sprintf("%s/dataset/%s/%s/metadata.json", output_final, point_name, basename(s2_id))

    dir.create(sprintf("%s/input", output_final_folder), showWarnings = FALSE, recursive = TRUE)
    dir.create(sprintf("%s/target", output_final_folder), showWarnings = FALSE, recursive = TRUE)
    dir.create(sprintf("%s/thumbnails", output_final_folder), showWarnings = FALSE, recursive = TRUE)

    # 4.9 Create cloud-segmentation.json (main file in iris software)
    metadata_main <- sprintf("%s/dataset/%s/cloud_segmentation_%s.json", output_final, point_name, point_name)
    if (counter == 0) {
      ee_create_cloudseg(path = metadata_main)
    }

    # 4.10 Users now have userID:password (lab:lab)
    if (counter == 0) {
       # lab_lab_user(path = dirname(metadata_main), point_name = point_name)
    }

    # 4.11 Generate a script to upload results to Google Drive database
    if (counter == 0) {
      generate_script(dirname(metadata_main))
    }

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
      IPL_multitemporal_cloud_logical = IPL_multitemporal_cloud_logical,
      output_final_folder = output_final_folder
    )

    # 4.16 Create metadata.json for each file
    ee_create_metadata(
      id = basename(s2_id),
      point = c(jsonfile_r$y, jsonfile_r$x),
      path = metadata_spec
    )

    # 4.17 Create a thumbnail
    ee_generate_thumbnail(
      s2_id = s2_id,
      final_stack =  final_stack,
      crs_kernel =  crs_kernel,
      output_final_folder = output_final_folder
    )
    counter <- counter + 1
  }

  row_position <- gsub("point_", "", point_name) %>% as.numeric()
  labelers_names <- c("Jhomira", "Fernando", "Eduardo")
  # point_metadata
  df_final <- data_frame(
    id = sprintf("%s_%02d", point_name, 1:5),
    labeler = sp_db[row_position,]$labeler,
    type = sp_db[row_position,]$label,
    difficulty = NA,
    sen2_id = s2_ids,
    sen2_date = s2_dates,
    sen1_id = s1_ids,
    s1_date = s1_dates,
    sen1_area = s1_area_per,
    land_use = land_use,
    elevation = elevation,
    shadow_dir = shadow_dir,
    split = NA,
    state = FALSE,
    evaluation_I = sp_db[row_position,]$validator,
    evaluation_II = labelers_names[!(labelers_names %in% c(sp_db[row_position,]$validator, sp_db[row_position,]$labeler))],
    evaluation_Expert = FALSE
  )

  write_csv(
    x = df_final,
    file = sprintf("%s/%s_metadata.csv", dirname(metadata_main), point_name)
  )

  # 5. Save geometry
  roi <- extent(final_stack[[1]]) %>%
    st_bbox() %>%
    st_as_sfc()
  st_crs(roi) <- crs_kernel
  write_sf(roi, sprintf("%s/%s.gpkg", dirname(output_final_folder), point_name))

  # Return is everything is fine :)
  invisible(TRUE)
}

#' Function to create metadata (i.e. metadata_1500.json)
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


#' Get the respective Sentinel-1 image for a specific Sentinel-2
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


#' Merge Sentinel-2, Sentinel-1 and cloud label (sen2cor, s2cloudness) predictions
#' @noRd
ee_merge_s2_full <- function(s2_id, s1_id, s2_date) {
  # 1. Define all the images that will be inside input.tif (ImageCollection)
  ## S1
  s1_grd <- ee$Image(s1_id)$unmask(-99, sameFootprint = F)
  ## S2-level 2A (sen2cor)
  s2_2a <- ee$Image(sprintf("COPERNICUS/S2_SR/%s", basename(s2_id)))
  ## S2-level 1C
  s2_1c <- ee$Image(sprintf("COPERNICUS/S2/%s", basename(s2_id)))
  ### Estimate CDI
  s2_cdi <- ee$Algorithms$Sentinel2$CDI(s2_1c) %>%
    ee$Image$unmask(-99, sameFootprint = F)
  ## DEM data (MERIT)
  extra_dem <- cloudsen12_dem()
  ## LandUSE data (copernicus)
  extra_LC <- cloudsen12_lc()

  # 2. Define all the images that will be inside target.tif (ImageCollection)
  ## Cloud mask according to Zupanc et al. 2019
  ## https://medium.com/sentinel-hub/sentinel-hub-cloud-detector-s2cloudless-a67d263d3025
  s2_cloud <- ee$Image(sprintf("COPERNICUS/S2_CLOUD_PROBABILITY/%s", basename(s2_id)))
  boxcar1 <- ee$Kernel$square(radius = 4, units = 'pixels')
  boxcar2 <- ee$Kernel$square(radius = 2, units = 'pixels')
  s2cloudness_prob <- s2_cloud$rename("cmask_s2cloudness")
  s2cloudness_prob_reclass <- s2_cloud %>%
    ee$Image$convolve(boxcar1) %>%
    ee$Image$focal_max(kernel = boxcar2, iterations = 1) %>%
    ee$Image$gte(40) %>%
    ee$Image$rename("cmask_s2cloudness_reclass")

  # 3. Reclass sen2cor (SLC, sentinel2 level2A)
  ## SEN2COR:
  ## 1	ff0004	Saturated or defective
  ## 2	868686	Dark Area Pixels
  ## 3	774b0a	Cloud Shadows
  ## 4	10d22c	Vegetation
  ## 5	ffff52	Bare Soils
  ## 6	0000ff	Water
  ## 7	818181	Clouds Low Probability / Unclassified
  ## 8	c0c0c0	Clouds Medium Probability
  ## 9	f1f1f1	Clouds High Probability
  ## 10	bac5eb	Cirrus
  ## 11	52fff9	Snow / Ice
  ##
  ## CLOUDSEN12 Reclass
  ##
  ## 2, 4, 5, 6, 7, 11 -> clear (0)
  ## 8, 9 -> thick cloud (1)
  ## 10 -> CIRRUS thin cloud (2)
  ## 3 -> cloud shadows (3)
  ## 0, 1 -> no data (4)
  ## OBS: In Earth Engine "no data" values are masked out.
  s2_scl <- s2_2a$select("SCL")$rename("cmask_sen2cor")
  s2_scl_reclass <- s2_2a$select("SCL")$remap(
    c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
    c(4, 4, 0, 3, 0, 0, 0, 0, 1, 1, 2, 0)
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

#' Create IRIS main json :)
#' @noRd
ee_create_cloudseg <- function(path) {
  point_name <- gsub("\\.json$","",paste0(strsplit(basename(path), "_")[[1]][3:4],collapse = "_"))
  cseg_list <- list(
    name = point_name,
    authentication_required = TRUE,
    images = list(
      path = list(
        B1 = "{id}/input/B1.tif",
        B2 = "{id}/input/B2.tif",
        B3 = "{id}/input/B3.tif",
        B4 = "{id}/input/B4.tif",
        B5 = "{id}/input/B5.tif",
        B6 = "{id}/input/B6.tif",
        B7 = "{id}/input/B7.tif",
        B8 = "{id}/input/B8.tif",
        B8A = "{id}/input/B8A.tif",
        B9 = "{id}/input/B9.tif",
        B10 = "{id}/input/B10.tif",
        B10_threshold = "{id}/input/B10_threshold.tif",
        B11 = "{id}/input/B11.tif",
        B12 = "{id}/input/B12.tif",
        CDI = "{id}/input/CDI.tif",
        VV = "{id}/input/VV.tif",
        VH = "{id}/input/VH.tif",
        angle = "{id}/input/angle.tif",
        elevation = "{id}/input/elevation.tif",
        landuse = "{id}/input/landuse.tif",
        s2cloudness = "{id}/target/s2cloudness_prob.tif",
        s2cloudness_reclass = "{id}/target/s2cloudness_reclass.tif",
        sen2cor = "{id}/target/sen2cor_real.tif",
        sen2cor_reclass = "{id}/target/sen2cor_reclass.tif"
      ),
      shape = c(511,511),
      metadata = "{id}/metadata.json",
      thumbnails = "{id}/thumbnails/thumbnail.png"
    ),
    segmentation = list(
      path = "{id}/target/manual.png",
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
        data = "$B10.B1**0.8*5",
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
        data = c("$B1.B1", "$B11.B1", "$B12.B1")
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
      CIRRUS_THRESHOLD = list(
        description = "Cirrus threshold > 0.025",
        type = "image",
        data = "$B10_threshold.B1"
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


#' Funcion de pato para subir a drive no es mia :X <3
#' @noRd
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
cloudsen12_dem <- function() {
  merit_dem <- ee$Image("MERIT/Hydro/v1_0_1")$select("elv")
  cryosat2_dem <- ee$Image("CPOM/CryoSat2/ANTARCTICA_DEM")$select("elevation")
  cloudsen12_dem <- ee$Image$add(merit_dem$unmask(0), cryosat2_dem$unmask(0))
  cloudsen12_dem
}


#' Proba-V-C3 + ANTARCTICA (70 -> snow)
#' @noRd
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

#' Search metajson in roy folder :)
#' @noRd
search_metajson <- function(pattern) {
  # Search metadata
  drive_jsonfile <- drive_ls(
    path = as_id("1fBGAjZkjPEpPr0p7c-LtJmfbLq3s87RK"),
    q = sprintf("name contains '%s'", pattern)
  )
  # Metadata no found
  if (nrow(drive_jsonfile) == 0) {
    stop("metadata no found")
  }
  #
  jsonfile <- try(
    drive_download(
      file = drive_jsonfile,
      path = paste0(tempdir(),"/", drive_jsonfile$name),
      overwrite = TRUE
    )
  )
  stop_5 <- 0
  while((class(jsonfile)[1] == "try-error") | stop_5 == 5) {
    jsonfile <- try(
      drive_download(
        file = drive_jsonfile,
        path = paste0(tempdir(),"/", drive_jsonfile$name),
        overwrite = TRUE
      )
    )
    stop_5 = stop_5 + 1
  }
  paste0(tempdir(),"/", drive_jsonfile$name)
}


# assign_imgs <- function(pnumber) {
#   set.seed(100)
#   labeler_names <- c("Jhomira", "Luis", "Eduardo")
#   img_labeler <- sapply(1:pnumber, function(x) sample(3, 3, replace = FALSE)) %>% as.numeric() %>% as.factor
#   levels(img_labeler) <- labeler_names
#   img_labeler
# }

#' Search all the metadata available in Roy folder
#' @noRd
drive_metadata_json <- function() {
  drive <- sprintf("%s/drive_dataset.Rdata", tempdir())
  drive_jsonfile <- drive_ls(as_id("1fBGAjZkjPEpPr0p7c-LtJmfbLq3s87RK"))
  save(drive_jsonfile, file = drive)
  drive_jsonfile
}


# calibration images
sen2_id <- c(
  "20181214T090351_20181214T092521_T33KTU",
  "20190309T053649_20190309T054326_T43SET",
  "20190402T083559_20190402T085915_T34KBC",
  "20190403T133229_20190403T133227_T22KFG",
  "20190623T092039_20190623T092036_T36VVK",
  "20190701T051659_20190701T052112_T44SME",
  "20190729T025549_20190729T025903_T50RKU",
  "20190828T161839_20190828T162908_T16RFV",
  "20190902T025541_20190902T025544_T52WFV",
  "20200125T145719_20200125T150550_T18LYH",
  "20200405T104619_20200405T105015_T31TBG",
  "20200410T142731_20200410T143742_T19HCA",
  "20200619T181919_20200619T181958_T12UVU",
  "20200712T130251_20200712T130252_T24LVQ",
  "20200721T180921_20200721T181627_T12TWR"
)

# Evaluation points
cloudsen12_point_validation <- c(
  "metadata_1382.json", "metadata_0299.json", "metadata_0466.json",
  "metadata_0470.json", "metadata_0503.json", "metadata_0823.json",
  "metadata_0838.json", "metadata_0903.json", "metadata_0908.json",
  "metadata_0919.json", "metadata_0985.json", "metadata_0995.json",
  "metadata_1004.json", "metadata_1031.json", "metadata_0183.json"
)

# detect_points <- function(x) {
#   fjson <- jsonlite::read_json(x) %>% names()
#   sum(sen2_id %in% fjson)
# }


#' Cloud probability (range of interest) to pick up images
#' @noRd
gen_rcloudpoints <- function(n) {
  # groups_n <- floor(c(0.05,0.3,0.3,0.3,0.05)*n)
  cloud_ppoints <- list()
  for (index in seq_len(n)) {
    groups_n <- rep(1, 5)
    # if (sum(groups_n) != n) {
    #   difff <- n - sum(groups_n)
    #   random_add <- sample(5,1)
    #   groups_n[random_add] <- groups_n[random_add] + difff
    # }
    cloud_ppoints[[index]] <- c(
      runif(n = groups_n[1], min = 0, max = 0),  # clear
      runif(n = groups_n[2], min = 5, max = 25), # almost clear
      runif(n = groups_n[3], min = 25, max = 45), # low-cloudy
      runif(n = groups_n[4], min = 45, max = 65), # mid-cloudy
      runif(n = groups_n[5], min = 65, max = 100)) # cloudy
  }
  cloud_ppoints %>%
    as.data.frame() %>%
    `colnames<-`(NULL) %>%
    t %>%
    as.data.frame()
}

# # # -------
# # # CloudSEN12: https://code.earthengine.google.com/d7a7751a50f18cefc66f18082a87eed3
# #
# # # ALGORITHM SETTINGS
# # cloudThresh <- 0.2 # Ranges from 0-1.Lower value will mask more pixels out. Generally 0.1-0.3 works well with 0.2 being used most commonly
# # cloudHeights <- ee$List$sequence(200,10000,250) # Height of clouds to use to project cloud shadows
# # irSumThresh <- 0.3 # Sum of IR bands to include as shadows within TDOM and the shadow shift method (lower number masks out less)
# # ndviThresh <- -0.1
# # dilatePixels <- 2 # Pixels to dilate around clouds
# # contractPixels <- 1 # Pixels to reduce cloud mask and dark shadows by to reduce inclusion of single-pixel comission errors
# # erodePixels <- 1.5
# # dilationPixels <- 3
# # cloudFreeKeepThresh <- 5
# # cloudMosaicThresh <- 50
#
#
# # calcCloudStats: Calculates a mask for clouds in the image.
# #        input: im - Image from image collection with a valid mask layer
# #        output: original image with added stats.
# #                - CLOUDY_PERCENTAGE: The percentage of the image area affected by clouds
# #                - ROI_COVERAGE_PERCENT: The percentage of the ROI region this particular image covers
# #                - CLOUDY_PERCENTAGE_ROI: The percentage of the original ROI which is affected by the clouds in this image
# #                - cloudScore: A per pixel score of cloudiness
# calcCloudStats <- function(img) {
#   imgPoly <- ee$Algorithms$GeometryConstructors$Polygon(
#     ee$Geometry(img$get("system:footprint"))$coordinates()
#   )
#   roi <- ee$Geometry(img$get("ROI"))
#   intersection <- roi$intersection(imgPoly, ee$ErrorMargin(0.5))
#   cloudMask <- img$select("cloudScore")$gt(cloudThresh)$clip(roi)$rename("cloudMask")
#   cloudAreaImg <- cloudMask$multiply(ee$Image$pixelArea())
#   stats <- cloudAreaImg$reduceRegion(
#     reducer = ee$Reducer$sum(),
#     geometry = roi,
#     scale = 10,
#     maxPixels = 1e12
#   )
#
#   cloudPercent <- ee$Number(stats$get("cloudMask"))$divide(imgPoly$area())$multiply(100)
#   coveragePercent <- ee$Number(intersection$area())$divide(roi$area())$multiply(100)
#   cloudPercentROI <- ee$Number(stats$get("cloudMask"))$divide(roi$area())$multiply(100)
#
#
#   img <- img$set("CLOUDY_PERCENTAGE", cloudPercent)
#   img <- img$set("ROI_COVERAGE_PERCENT", coveragePercent)
#   img <- img$set("CLOUDY_PERCENTAGE_ROI", cloudPercentROI)
#   img
# }
#
# rescale <- function(img, exp, thresholds) {
#   img$expression(exp, list(img = img))$
#     subtract(thresholds[1])$
#     divide(thresholds[2] - thresholds[1])
# }
#
#
# computeQualityScore <- function(img) {
#   score <- img$select("cloudScore")$max(img$select("shadowScore"))
#   score <- score$reproject("EPSG:4326", NULL, 20)$reduceNeighborhood(
#     reducer = ee$Reducer$mean(),
#     kernel = ee$Kernel$square(5)
#   )
#   score <- score$multiply(-1)
#   img$addBands(score$rename("cloudShadowScore"))
# }
#
#
# #**
# # Implementation of Basic cloud shadow shift
# #
# # Author: Gennadii Donchyts
# # License: Apache 2.0
# #
# projectShadows <- function(image) {
#   meanAzimuth <- image$get("MEAN_SOLAR_AZIMUTH_ANGLE")
#   meanZenith <- image$get("MEAN_SOLAR_ZENITH_ANGLE")
#   cloudMask <- image$select("cmask_s2cloudness_reclass")
#
#   #Find dark pixels
#   darkPixelsImg <- image$select(c("B8", "B11", "B12"))$
#     divide(10000)$
#     reduce(ee$Reducer$sum())
#
#   ndvi <- image$normalizedDifference(c("B8", "B4"))
#   waterMask <- ndvi$lt(ndviThresh)
#   darkPixels <- darkPixelsImg$lt(irSumThresh)
#
#   # Get the mask of pixels which might be shadows excluding water
#   darkPixelMask <- darkPixels$And(waterMask$Not())
#   darkPixelMask <- darkPixelMask$And(cloudMask$Not())
#
#   #Find where cloud shadows should be based on solar geometry
#   #Convert to radians
#   azR <- ee$Number(meanAzimuth)$add(180)$multiply(base::pi)$divide(180.0)
#   zenR <- ee$Number(meanZenith)$multiply(base::pi)$divide(180.0)
#
#   # Find the shadows
#   # cloudHeights <- ee$List$sequence(200,10000,1000) # Height of clouds to use to project cloud shadows
#   cloudHeight <- 1000
#   shadows <- cloudHeights$map(
#     ee_utils_pyfunc(
#       function(cloudHeight){
#         cloudHeight <- ee$Number(cloudHeight)
#         shadowCastedDistance <- zenR$tan()$multiply(cloudHeight) # Distance shadow is cast
#         x <- azR$sin()$multiply(shadowCastedDistance)$multiply(-1) #.divide(nominalScale); #X distance of shadow
#         y <- azR$cos()$multiply(shadowCastedDistance)$multiply(-1) #Y distance of shadow
#         image$select("cmask_s2cloudness")$displace(ee$Image$constant(x)$addBands(ee$Image$constant(y)))
#       }
#     )
#   )
#
#   shadowMasks <- ee$ImageCollection$fromImages(shadows)
#   shadowMask <- shadowMasks$mean()
#
#
#   # Create shadow mask
#   shadowMask <- dilatedErossion(shadowMask$multiply(darkPixelMask))
#   shadowScore <- shadowMask$reduceNeighborhood(
#     reducer = ee$Reducer$max(),
#     kernel = ee$Kernel$square(1)
#   )
#   image <- image$addBands(shadowScore$rename("shadowScore"))
#   image
# }
#
# dilatedErossion <- function(score) {
#   score$reproject('EPSG:4326', NULL, 20)$
#     focal_min(radius = erodePixels, iterations = 3)$
#     focal_max(radius = dilationPixels, iterations = 3)$
#     reproject("EPSG:4326", NULL, 20)
# }
#
# computeS2CloudScore <- function(img) {
#   toa <- img$select(c('B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12'))$
#     divide(10000)
#   toa <- toa$addBands(img$select("QA60"))
#
#   # ['QA60', 'B1','B2',    'B3',    'B4',   'B5','B6','B7', 'B8','  B8A', 'B9', 'B10', 'B11','B12']
#   # ['QA60','cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2']);
#   # Compute several indicators of cloudyness and take the minimum of them.
#   score <- ee$Image(1)
#
#   # Clouds are reasonably bright in the blue and cirrus bands.
#   score <- score$min(rescale(toa, 'img.B2', c(0.1, 0.5)))
#   score <- score$min(rescale(toa, 'img.B1', c(0.1, 0.3)))
#   score <- score$min(rescale(toa, 'img.B1 + img.B10', c(0.15, 0.2)))
#
#   # Clouds are reasonably bright in all visible bands.
#   score <- score$min(rescale(toa, 'img.B4 + img.B3 + img.B2', c(0.2, 0.8)))
#
#   # Clouds are moist
#   ndmi <- img$normalizedDifference(c("B8", "B11"))
#   score <- score$min(rescale(ndmi, "img", c(-0.1, 0.1)))
#
#   # However, clouds are not snow.
#   ndsi <- img$normalizedDifference(c("B3", "B11"))
#   score <- score$min(rescale(ndsi, 'img', c(0.8, 0.6)))
#
#   # Clip the lower end of the score
#   score <- score$max(ee$Image(0.001))
#
#   # Remove small regions and clip the upper bound
#   dilated <- dilatedErossion(score)$min(ee$Image(1.0))
#
#   # score = score.multiply(dilated)
#   score <- score$reduceNeighborhood(
#     reducer = ee$Reducer$mean(),
#     kernel = ee$Kernel$square(5)
#   )
#
#   img$addBands(score$rename("cloudScore"))
# }


#' IPL multitemporal cloud detection algorithm
#' Issue: This algorithm required images with a CC < 5%. This condition not always
#' can be fulfilled. ee_upl_cloud_logical test if there are enough images
#' to run the algorithm.
#' @param sen2id The ID of the sentinel2 image to test
#' @param roi Region of interest
ee_upl_cloud_logical <- function(sen2id, roi) {
  image_collection_name <- 'COPERNICUS/S2/'
  image_wrap <- ee_cloud$image_wrapper$S2L1CImage(
    collection = image_collection_name,
    index = sen2id
  )
  cloud_img_cc <- ee_cloud$multitemporal_cloud_masking$CloudClusterScore(
    img = image_wrap,
    region_of_interest = roi,
    method_pred="persistence"
  )[[1]]

  ee_results_cloud <- try(expr = cloud_img_cc$getInfo(), silent = TRUE)
  if (class(ee_results_cloud) == "try-error") {
    FALSE
  } else {
    TRUE
  }
}

#' IPL multitemporal cloud detection algorithm
#' @param sen2id The ID of the sentinel2 image to test
#' @param roi Region of interest
#' @noRd
ee_upl_cloud <- function(sen2id, roi) {
  image_collection_name <- 'COPERNICUS/S2/'
  image_wrap <- ee_cloud$image_wrapper$S2L1CImage(
    collection = image_collection_name,
    index = sen2id
  )
  cloud_img_cc <- ee_cloud$multitemporal_cloud_masking$CloudClusterScore(
    img = image_wrap,
    region_of_interest = roi,
    method_pred="persistence"
  )[[1]]
  cloud_img_cc %>%
    ee$Image$rename("IPL_cloud_mask") %>%
    ee$Image$remap(
      c(0, 1, 2),
      c(0, 3, 1)
    )
}

#' Estimate the shadow direction considering the SOLAR ZENITH and SOLAR AZIMUTH.
#' @param image A Sentinel2 image
#' @noRd
shadow_direction <- function(image) {
  meanAzimuth <- image$get("MEAN_SOLAR_AZIMUTH_ANGLE")
  meanZenith <- image$get("MEAN_SOLAR_ZENITH_ANGLE")
  azR <- ee$Number(meanAzimuth)$add(180)$multiply(base::pi)$divide(180.0)
  zenR <- ee$Number(meanZenith)$multiply(base::pi)$divide(180.0)
  shadowCastedDistance <- zenR$tan()
  x <- azR$sin()$multiply(shadowCastedDistance)
  y <- azR$cos()$multiply(shadowCastedDistance)
  ee$Number$atan2(x, y)$multiply(180)$divide(base::pi) %>%
    ee$Image() %>%
    ee$Image$rename("cloudshadow_direction")
}


#' Add 360 if the shadow direction is negative
#' @param raster RasterStack
fix_shadow_direction <- function(raster) {
  if (mean(getValues(raster)) < 0) {
    raster <- raster + 360
  } else {
    raster
  }
}

#' Save Cloud mask available in Earth Engine
#' @noRd
create_target_raster <- function(final_stack, IPL_multitemporal_cloud_logical, output_final_folder) {
  if (!IPL_multitemporal_cloud_logical) {
    bandnames <- c("s2cloudness_prob", "s2cloudness_reclass", "sen2cor_real", "sen2cor_reclass")
    benchmarch_data <- stack(final_stack[[c(18:21)]])
    target_spec <- sprintf("%s/target/%s.tif", output_final_folder, bandnames)
    lapply(1:4, function(x) writeRaster(benchmarch_data[[x]], target_spec[x], overwrite = TRUE))
  } else {
    IPL_cloudmask_reclass <- final_stack[[25]]
    benchmarch_data <- stack(final_stack[[c(18:21)]], IPL_cloudmask_reclass)
    bandnames <- c("s2cloudness_prob", "s2cloudness_reclass", "sen2cor_real", "sen2cor_reclass", "IPL_cloudmask_reclass")
    target_spec <- sprintf("%s/target/%s.tif", output_final_folder, bandnames)
    lapply(1:5, function(x) writeRaster(benchmarch_data[[x]], target_spec[x], overwrite = TRUE))
  }
}

#' Create and save a thumbnail
#' @noRd
ee_generate_thumbnail <- function(s2_id, final_stack, crs_kernel, output_final_folder) {
  # Get area stk - polygon
  extent_stk <- extent(final_stack) %>%
    as("SpatialPolygons") %>%
    st_as_sfc()
  st_crs(extent_stk) <- crs_kernel

  # Double - polygon
  thumbnail_buffer <- st_buffer(
    x = st_centroid(extent_stk),
    dist = 10 * 510,
    endCapStyle = "SQUARE",
    joinStyle = "BEVEL"
  )

  # Download thumbnail
  s2_img432 <- ee_as_raster(
    image = ee$Image(s2_id)$select(c("B4","B3","B2"))$unitScale(0, 10000)$multiply(255)$toByte(),
    scale = 10,
    region = sf_as_ee(thumbnail_buffer),
    quiet = TRUE
  )

  # Save thumbnail tif
  writeRaster(s2_img432, sprintf("%s/thumbnails/thumbnail.tif", output_final_folder), overwrite = TRUE)

  # Save thumbnail png
  png(sprintf("%s/thumbnails/thumbnail.png", output_final_folder), 1021, 1021)
  plotRGB(s2_img432, stretch = "lin")
  plot(extent_stk, add=TRUE, border = "red", lwd = 3)
  on.exit(dev.off())
}


#' Create lab:lab user
lab_lab_user <- function(path, point_name) {
  dir.create(
    sprintf("%s/%s.iris", path, point_name),
    recursive = TRUE,
    showWarnings = FALSE
  )
  download.file(
    "https://drive.google.com/u/1/uc?id=1BZB3ZIcRnpB0wIEFJCXE5XsH7mh2gI7D&export=download",
    sprintf("%s/%s.iris/iris.db", path, point_name)
  )
}



#' Dilate high-quality results using a 3x3 kernel
#' @noRd
dilate_raster <- function(label_raster) {
  # Labels according to cloudsen12
  # 0 -> clear
  # 1 -> thick cloud
  # 2 -> thin cloud
  # 3 -> shadow

  # Re-order raster value
  clear <- (label_raster == 0) * 0
  thick_cloud <- (label_raster == 1) * 3
  thin_cloud <- (label_raster == 2) * 1
  cloud_shadow <- (label_raster == 3) * 2

  order_new_r <- clear + thick_cloud + thin_cloud + cloud_shadow
  order_new_r_values <- as.matrix(order_new_r)

  # 3x3 kernel
  (k <- matrix(1,nrow=3, ncol=3))

  # dilate raster R
  rdilate <- dilate(order_new_r_values, k)
  order_new_r[] <- rdilate


  # Re-order raster value
  clear <- (order_new_r == 0) * 0
  thick_cloud <- (order_new_r == 3) * 1
  thin_cloud <- (order_new_r == 1) * 2
  cloud_shadow <- (order_new_r == 2) * 3

  # final_raster
  final_r <- clear + thick_cloud + thin_cloud + cloud_shadow
  final_r
}


#' Generate a preview in a specific image
#' @noRd
generate_spplot <-function(rgb, label_raster, output = "comparison.svg") {
  # Convert RGB to use with spplot
  lout <- rgb2spLayout(rgb)

  # raster label remove clear
  label_raster[label_raster == 0] = NA

  # all is NA?
  rat <- nrow(levels(as.factor(label_raster))[[1]])

  # create viz folder
  dir.create(path = "viz", showWarnings = FALSE)

  if (rat > 0) {
    # label spplot
    ## FFFFFF <- Clear (0)
    ## FFFF00 <- Thick cloud (1)
    ## 00FF00 <- Thin cloud (2)
    ## FF0000 <- Shadow (3)
    spplot1 <- spplot(
      obj = label_raster,
      sp.layout = list(
        lout
      ),
      col.regions = c("#FFFF00", "#00FF00", "#FF0000"),
      at = c(0,1,2,3,4),
      par.settings = list(axis.line = list(col = "transparent")),
      colorkey=FALSE
    )
  }

  fake_r <- label_raster
  fake_r[fake_r >= 0] = NA
  fake_r[1,1] = 0

  # RGB spplot
  spplot2 <- spplot(
    obj = fake_r,
    sp.layout = lout,
    cex = 0,
    par.settings = list(axis.line = list(col = "transparent")),
    colorkey=FALSE,
    alpha.regions = 0
  )
  if (rat > 0) {
    svg(output, width = 7*2, height = 7*2)
    grid.arrange(spplot2, spplot1, nrow = 1)
  } else {
    svg(output, width = 7*2, height = 7*2)
    grid.arrange(spplot2, spplot2, nrow=1)
  }
  on.exit(dev.off())
}

#' Hex color (manual.png) to raster
#' @noRd
number_to_hex <- function(stack_color, raster_ref) {
  raster_ref_f <- raster_ref
  stk_r <- stack(stack_color)

  #FFFFFF <- Clear (0)
  #FFFF00 <- Thick cloud (1)
  #00FF00 <- Thin cloud (2)
  #FF0000 <- Shadow (3)
  color_values <- rgb(stk_r[],maxColorValue = 255)
  color_values[color_values == "#FFFFFF"] <- 0
  color_values[color_values == "#FFFF00"] <- 1
  color_values[color_values == "#00FF00"] <- 2
  color_values[color_values == "#FF0000"] <- 3
  color_values <- as.numeric(color_values)

  raster_ref_f[] <- color_values
  raster_ref_f
}

#' Generate previews (in the 5 images)
#' @noRd
upload_results <- function(download_csaybar_user = FALSE) {

  if (download_csaybar_user) {
    dir_p  <- paste0(tempdir(), "/cd26ed5dc626f11802a652e81d02762e_s1078735@stud.sbg.ac.at")
    #download privileges
    download.file(
      url = "https://drive.google.com/uc?id=1dg02Ue6rkFa7xn7TU57bWu_gQxuxpGdY&export=download",
      destfile = dir_p
    )

    drive_auth("s1078735@stud.sbg.ac.at", token = dir_p)
  }

  # List files (Get SENTINEL_2 ID)
  point_name <- gsub("\\.gpkg$", "",list.files(pattern = "\\.gpkg$"))

  # Points Google Drive ID
  files_points_general <- googledrive::drive_ls(
    path = as_id("1BeVp0i-dGSuBqCQgdGZVDj4qzX1ms7L6"),
    q = sprintf("name contains '%s'", point_name)
  )

  # files Google Drive ID
  files_points <- googledrive::drive_ls(
    path = as_id(files_points_general$id)
  )
  folder_id <- files_points_general$id

  # List files (Get SENTINEL_2 ID)
  directories <- list.files()
  directories_sen2ids <- directories[grepl("^[0-9]", directories)]

  # Get the point number
  point_name <- gsub("\\.gpkg$", "", directories[grepl("gpkg", directories)])

  # Reference raster (base_raster)
  reference_raster <- raster(paste0(directories_sen2ids[1],"/input/B1.tif"))

  for (directory in directories_sen2ids) {
    message("Uploading: ",directory)
    # From png to raster
    high_labeling_r <- number_to_hex(
      stack_color = paste0(directory,"/target/manual.png"),
      raster_ref = reference_raster
    )

    # Dilate raster 3x3
    high_labeling_r_dilate <- high_labeling_r
    tmp_r <- paste0(tempdir(), "/manual.tif")
    writeRaster(x = high_labeling_r_dilate, filename = tmp_r, overwrite = TRUE)

    # Upload Raster to Google Drive
    folder_gd_img <- files_points[files_points$name %in% directory,]$id
    folder_gd_img_target <- googledrive::drive_ls(
      path = as_id(folder_gd_img)
    )
    target_ID <- folder_gd_img_target[folder_gd_img_target$name == "target", ]$id
    drive_upload(
      media = tmp_r,
      overwrite = TRUE,
      path = as_id(target_ID),
      verbose = FALSE
    )
  }

  # point.iris ID folder
  point_iris <- files_points[grepl("iris", files_points$name), ]$id

  # Create zip file
  zip::zip(
    zipfile = paste0(point_name, ".zip"),
    files = list.files(
      path = paste0(point_name, ".iris/segmentation"),
      recursive = TRUE,
      full.names = TRUE
    )
  )

  # Upload zip file - point
  drive_upload(
    media = paste0(point_name, ".zip"),
    overwrite = TRUE,
    path = as_id(folder_id),
    verbose = FALSE
  )

  # Upload metadata - metadata
  drive_upload(
    media = paste0(point_name, "_metadata.csv"),
    overwrite = TRUE,
    path = as_id(folder_id),
    verbose = FALSE
  )

  # Upload  - VIZ zip
  drive_upload(
    media = paste0("viz/",point_name, "_viz.zip"),
    overwrite = TRUE,
    path = as_id(folder_id),
    verbose = FALSE
  )
}

#' generate plots
#' @noRd
generate_preview <- function() {
  # Point name
  point_name <- gsub("\\.gpkg$", "",list.files(pattern = "\\.gpkg$"))

  # List files (Get SENTINEL_2 ID)
  directories <- list.files()
  directories_sen2ids <- directories[grepl("^[0-9]", directories)]

  # Create Viz folder
  dir.create(path = "viz", showWarnings = FALSE)

  # Get the point number
  for (directory in directories_sen2ids) {
    # RGB -> read STACK RasterStack
    rgb <- stack(sprintf(paste0(directory, "/input/B%01d.tif"), 4:2))

    # from .png to rasterlayer
    high_labeling_r <- number_to_hex(
      stack_color = paste0(directory,"/target/manual.png"),
      raster_ref = rgb[[1]]
    ) #%>% dilate_raster()

    # Generate .svg
    generate_spplot(
      rgb = rgb,
      label_raster = high_labeling_r,
      output = sprintf("viz/%s.svg", directory)
    )
  }

  # svg files
  svg_files <- list.files("viz", "\\.svg$",full.names = T)

  # Create zip file
  zip::zip(
    zipfile = paste0("viz/",point_name, "_viz.zip"),
    files = list.files(
      path = "viz",
      recursive = TRUE,
      full.names = TRUE
    )
  )
}


#' Download labelers
#' @noRd
download_labels <- function(point) {
  # 1. List files (Get SENTINEL_2 ID)
  point_name <- point

  # 2. Create dir
  dir.create(
    path = sprintf("../%s/%s.iris/segmentation", point, point),
    showWarnings = FALSE,
    recursive = TRUE
  )


  # 2. Points Google Drive ID
  files_points_general <- googledrive::drive_ls(
    path = as_id("1BeVp0i-dGSuBqCQgdGZVDj4qzX1ms7L6"),
    q = sprintf("name contains '%s'", point_name)
  )

  # 3. Files inside the folder in GoogleDrive
  files_points <- googledrive::drive_ls(
    path = as_id(files_points_general$id)
  )

  zip_files <- files_points[files_points$name == paste0(point_name, ".zip"),]$id

  drive_download(
    file = as_id(zip_files),
    path = paste0(point_name, ".zip"),
    overwrite = TRUE
  )

  # 4. Unzip file
  zip::unzip(
    zipfile = paste0(point_name, ".zip"),
    exdir = sprintf("../%s", point)
  )
}


# Download viz
download_viz <- function(point) {
  # 1. List files (Get SENTINEL_2 ID)
  point_name <- point
  point_name_viz <- paste0(point_name, "_viz.zip")

  # Create dir
  dir.create(
    path = paste0("../", point,"/viz"),
    showWarnings = FALSE,
    recursive = TRUE
  )

  # Points Google Drive ID
  files_points_general <- googledrive::drive_ls(
    path = as_id("1BeVp0i-dGSuBqCQgdGZVDj4qzX1ms7L6"),
    q = sprintf("name contains '%s'", point_name)
  )

  # files Google Drive ID
  files_points <- googledrive::drive_ls(
    path = as_id(files_points_general$id)
  )

  # Download zip
  googledrive::drive_download(
    file = as_id(files_points[files_points$name == point_name_viz,]$id),
    path = paste0("../", point,"/viz/", point_name_viz),
    overwrite = TRUE
  )
  #  UNzip
  zip::unzip(
    zipfile = paste0("../", point,"/viz/", point_name_viz),
    exdir = paste0("../", point)
  )
}

#' Upload preview
#' @noRd
upload_preview <- function(download_csaybar_user = FALSE) {
  if (download_csaybar_user) {
    # csaybar user
    dir_p  <- paste0(tempdir(), "/cd26ed5dc626f11802a652e81d02762e_s1078735@stud.sbg.ac.at")

    # download privileges
    download.file(
      url = "https://drive.google.com/uc?id=1dg02Ue6rkFa7xn7TU57bWu_gQxuxpGdY&export=download",
      destfile = dir_p
    )

    # Auth Google Drive
    drive_auth("s1078735@stud.sbg.ac.at", token = dir_p)
  }

  # List files (Get SENTINEL_2 ID)
  point_name <- gsub("\\.gpkg$", "",list.files(pattern = "\\.gpkg$"))

  # Points Google Drive ID
  files_points_general <- googledrive::drive_ls(
    path = as_id("1BeVp0i-dGSuBqCQgdGZVDj4qzX1ms7L6"),
    q = sprintf("name contains '%s'", point_name)
  )

  # files Google Drive ID
  files_points <- googledrive::drive_ls(
    path = as_id(files_points_general$id)
  )
  folder_id <- files_points_general$id

  # List files (Get SENTINEL_2 ID)
  svg_files <- list.files("viz", "\\.svg$",full.names = T)

  # Create zip file
  zip::zip(
    zipfile = paste0("viz/",point_name, "_viz.zip"),
    files = list.files(
      path = "viz",
      recursive = TRUE,
      full.names = TRUE
    )
  )

  # Upload  - VIZ zip
  drive_upload(
    media = paste0("viz/",point_name, "_viz.zip"),
    overwrite = TRUE,
    path = as_id(folder_id),
    verbose = FALSE
  )
}

#' Download thumbnails
#' @noRd
download_thumbnails <- function(point) {
  # 1. List files (Get SENTINEL_2 ID)
  point_name <- point

  # Create dir
  dir.create(
    path = paste0("../", point,"/thumbnails"),
    showWarnings = FALSE,
    recursive = TRUE
  )

  # Points Google Drive ID
  files_points_general <- googledrive::drive_ls(
    path = as_id("1BeVp0i-dGSuBqCQgdGZVDj4qzX1ms7L6"),
    q = sprintf("name contains '%s'", point_name)
  )

  # files Google Drive ID
  files_points <- googledrive::drive_ls(
    path = as_id(files_points_general$id),
    recursive = TRUE
  )

  # Download thumbnails.png
  thumbnails_to_download <- files_points[files_points$name %in% "thumbnail.png",]
  thumbnails_to_download$name <- sprintf("thumbnails_%02d.png", 1:5)

  for (index in 1:5) {
    googledrive::drive_download(
      file = thumbnails_to_download[index,],
      path = paste0("../", point,"/thumbnails/", thumbnails_to_download$name)[index],
      overwrite = TRUE
    )
  }
}

#' Generate R script
#' @noRd
generate_script <- function(path) {
  fileConn <- file(paste0(path, "/Run.R"))
  writeLines(
    text = c(
      "library(googledrive)",
      "library(tidyverse)",
      "library(gridExtra)",
      "library(raster)",
      "library(mmand)",
      "library(Orcs)",
      "library(zip)",
      "",
      "source(\"https://gist.githubusercontent.com/csaybar/daa1a877f3d1703b61846603e986b14c/raw/735282528b8088d58af49c4f295d54ae7b7daa5c/demo.R\")",
      "",
      "# Generar un .svg para cada imagen.",
      "generate_preview()",
      "",
      "# Subir todos los resultados de Google Drive",
      "httr::set_config(httr::config( ssl_verifypeer = 0L))",
      "httr::set_config(httr::config(http_version = 0))",
      "upload_results()",
      "",
      "# Subir solo los svgs",
      "upload_preview()",
      "",
      "",
      "# Funciones auxiliares para descargar",
      "# download_viz(point = \"point_1382\")",
      "# download_thumbnails(point = \"point_1382\")",
      "# download_labels(point = \"point_1382\")"
    ),
    con =  fileConn
  )
  close(fileConn)
}


# Mode in R
Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

#' Dilate high-quality results using a 3x3 kernel
#' @noRd
dilate_raster <- function(label_raster) {
  # Labels according to cloudsen12
  # 0 -> clear
  # 1 -> thick cloud
  # 2 -> thin cloud
  # 3 -> shadow

  # Re-order raster value
  clear <- (label_raster == 0) * 0
  thick_cloud <- (label_raster == 1) * 3
  thin_cloud <- (label_raster == 2) * 1
  cloud_shadow <- (label_raster == 3) * 2

  order_new_r <- clear + thick_cloud + thin_cloud + cloud_shadow
  order_new_r_values <- as.matrix(order_new_r)

  # 3x3 kernel
  (k <- matrix(1,nrow=3, ncol=3))

  # dilate raster R
  rdilate <- dilate(order_new_r_values, k)
  order_new_r[] <- rdilate


  # Re-order raster value
  clear <- (order_new_r == 0) * 0
  thick_cloud <- (order_new_r == 3) * 1
  thin_cloud <- (order_new_r == 1) * 2
  cloud_shadow <- (order_new_r == 2) * 3

  # final_raster
  final_r <- clear + thick_cloud + thin_cloud + cloud_shadow
  final_r
}

#' Generate a preview in a specific image
#' @noRd
generate_spplot <-function(rgb, label_raster, output = "comparison.svg") {
  # Convert RGB to use with spplot
  lout <- rgb2spLayout(rgb)

  # raster label remove clear
  label_raster[label_raster == 0] = NA

  # all is NA?
  rat <- nrow(levels(as.factor(label_raster))[[1]])

  # create viz folder
  dir.create(path = "viz", showWarnings = FALSE)

  if (rat > 0) {
    # label spplot
    ## FFFFFF <- Clear (0)
    ## FFFF00 <- Thick cloud (1)
    ## 00FF00 <- Thin cloud (2)
    ## FF0000 <- Shadow (3)
    spplot1 <- spplot(
      obj = label_raster,
      sp.layout = list(
        lout
      ),
      col.regions = c("#FFFF00", "#00FF00", "#FF0000"),
      at = c(0,1,2,3,4),
      par.settings = list(axis.line = list(col = "transparent")),
      colorkey=FALSE
    )
  }

  fake_r <- label_raster
  fake_r[fake_r >= 0] = NA
  fake_r[1,1] = 0

  # RGB spplot
  spplot2 <- spplot(
    obj = fake_r,
    sp.layout = lout,
    cex = 0,
    par.settings = list(axis.line = list(col = "transparent")),
    colorkey=FALSE,
    alpha.regions = 0
  )
  if (rat > 0) {
    svg(output, width = 7*2, height = 7*2)
    grid.arrange(spplot2, spplot1, nrow = 1)
  } else {
    svg(output, width = 7*2, height = 7*2)
    grid.arrange(spplot2, spplot2, nrow=1)
  }
  on.exit(dev.off())
}

#' Hex color (manual.png) to raster
number_to_hex <- function(stack_color, raster_ref) {
  raster_ref_f <- raster_ref
  stk_r <- stack(stack_color)

  #FFFFFF <- Clear (0)
  #FFFF00 <- Thick cloud (1)
  #00FF00 <- Thin cloud (2)
  #FF0000 <- Shadow (3)
  color_values <- rgb(stk_r[],maxColorValue = 255)
  color_values[color_values == "#FFFFFF"] <- 0
  color_values[color_values == "#FFFF00"] <- 1
  color_values[color_values == "#00FF00"] <- 2
  color_values[color_values == "#FF0000"] <- 3
  color_values <- as.numeric(color_values)

  raster_ref_f[] <- color_values
  raster_ref_f
}

#' Generate previews (in the 5 images)
#' @noRd
upload_results <- function() {

  dir_p  <- paste0(tempdir(), "/cd26ed5dc626f11802a652e81d02762e_s1078735@stud.sbg.ac.at")
  #download privileges
  download.file(
    url = "https://drive.google.com/uc?id=1dg02Ue6rkFa7xn7TU57bWu_gQxuxpGdY&export=download",
    destfile = dir_p
  )

  drive_auth("s1078735@stud.sbg.ac.at", token = dir_p)

  # List files (Get SENTINEL_2 ID)
  point_name <- gsub("\\.gpkg$", "",list.files(pattern = "\\.gpkg$"))

  # Points Google Drive ID
  files_points_general <- googledrive::drive_ls(
    path = as_id("1BeVp0i-dGSuBqCQgdGZVDj4qzX1ms7L6"),
    q = sprintf("name contains '%s'", point_name)
  )

  # files Google Drive ID
  files_points <- googledrive::drive_ls(
    path = as_id(files_points_general$id)
  )
  folder_id <- files_points_general$id

  # List files (Get SENTINEL_2 ID)
  directories <- list.files()
  directories_sen2ids <- directories[!grepl("iris|viz", directories)]
  directories_sen2ids <- directories_sen2ids[dir.exists(directories_sen2ids)]

  # Get the point number
  point_name <- gsub("\\.gpkg$", "", directories[grepl("gpkg", directories)])

  # Reference raster (base_raster)
  reference_raster <- raster(paste0(directories_sen2ids[1],"/input/B1.tif"))

  for (directory in directories_sen2ids) {
    message("Uploading: ",directory)
    # From png to raster
    high_labeling_r <- number_to_hex(
      stack_color = paste0(directory,"/target/manual.png"),
      raster_ref = reference_raster
    )

    # Dilate raster 3x3
    high_labeling_r_dilate <- dilate_raster(high_labeling_r)
    tmp_r <- paste0(tempdir(), "/manual.tif")
    writeRaster(x = high_labeling_r_dilate, filename = tmp_r, overwrite = TRUE)

    # Upload Raster to Google Drive
    folder_gd_img <- files_points[files_points$name %in% directory,]$id
    folder_gd_img_target <- googledrive::drive_ls(
      path = as_id(folder_gd_img)
    )
    target_ID <- folder_gd_img_target[folder_gd_img_target$name == "target", ]$id
    drive_upload(
      media = tmp_r,
      overwrite = TRUE,
      path = as_id(target_ID),
      verbose = FALSE
    )
  }

  # point.iris ID folder
  point_iris <- files_points[grepl("iris", files_points$name), ]$id

  # Create zip file
  zip(
    zipfile = paste0(point_name, ".zip"),
    files = list.files(
      path = paste0(point_name, ".iris/segmentation"),
      recursive = TRUE,
      full.names = TRUE
    ),
    flags = "-q"
  )

  # Upload zip file
  drive_upload(
    media = paste0(point_name, ".zip"),
    overwrite = TRUE,
    path = as_id(folder_id),
    verbose = FALSE
  )

  # Upload metadata
  drive_upload(
    media = paste0(point_name, "_metadata.csv"),
    overwrite = TRUE,
    path = as_id(folder_id),
    verbose = FALSE
  )
}

#' generate plots
#' @noRd
generate_preview <- function() {
  # List files (Get SENTINEL_2 ID)
  directories <- list.files()
  directories_sen2ids <- directories[!grepl("iris", directories)]
  directories_sen2ids <- directories_sen2ids[dir.exists(directories_sen2ids)]
  directories_sen2ids <- directories_sen2ids[!directories_sen2ids == "viz"]

  # Create Viz folder
  dir.create(path = "viz", showWarnings = FALSE)

  # Get the point number
  for (directory in directories_sen2ids) {
    # RGB -> read STACK RasterStack
    rgb <- stack(sprintf(paste0(directory, "/input/B%01d.tif"), 4:2))

    # from .png to rasterlayer
    high_labeling_r <- number_to_hex(
      stack_color = paste0(directory,"/target/manual.png"),
      raster_ref = rgb[[1]]
    ) %>% dilate_raster()

    # Generate .svg
    generate_spplot(
      rgb = rgb,
      label_raster = high_labeling_r,
      output = sprintf("viz/%s.svg", directory)
    )
  }
}

# Download labelers
download_labels <- function() {
  # 1. List files (Get SENTINEL_2 ID)
  point_name <- gsub("\\.gpkg$", "",list.files(pattern = "\\.gpkg$"))

  # 2. Points Google Drive ID
  files_points_general <- googledrive::drive_ls(
    path = as_id("1BeVp0i-dGSuBqCQgdGZVDj4qzX1ms7L6"),
    q = sprintf("name contains '%s'", point_name)
  )

  # 3. Files inside the folder in GoogleDrive
  files_points <- googledrive::drive_ls(
    path = as_id(files_points_general$id)
  )
  zip_files <- files_points[files_points$name == paste0(point_name, ".zip"),]$id
  drive_download(
    file = as_id(zip_files),
    path = paste0(point_name, ".zip"),
    overwrite = TRUE
  )

  # 4. Unzip file
  unzip(paste0(point_name, ".zip"))
}


# Full download images
download_cloudSEN12_images <- function(points, local_cloudsen2_points, output) {
  # Create folder
  if (!dir.exists(output)) {
    dir.create(path = output, recursive = TRUE, showWarnings = TRUE)
  }

  for (index in points) {
    message("Downloading: Point_", index)
    metadata_json <- sprintf("metadata_%04d.json", index)
    jsonfile <- try(search_metajson(pattern = metadata_json))
    if (class(jsonfile) != "try-error") {
      results <- try(
        dataset_creator_chips(
          jsonfile = jsonfile,
          sp_db = local_cloudsen2_points,
          output_final = output
        )
      )
    }
  }
}

# Full download image thumbnails
select_dataset_thumbnail_creator_batch <- function(points,
                                                   local_cloudsen2_points,
                                                   n_images = "max",
                                                   kernel_size = c(255, 255),
                                                   data_range = c("2018-01-01", "2020-07-31"),
                                                   output = "results/") {
  for (index in points) {
    cloudsen2_row <- local_cloudsen2_points[index,]
    select_dataset_thumbnail_creator(
      cloudsen2_row = cloudsen2_row,
      n_images = n_images,
      kernel_size = kernel_size,
      data_range = data_range,
      output = output
    )
  }
}

# -------------------------------------------------------------------------
# MIGRATION FUNCTIONS -----------------------------------------------------
# -------------------------------------------------------------------------

landuse_types_list <- list(
  'Unknown' = 0,
  'Shrubs' = 20,
  'Herbaceous vegetation' = 30,
  'Cultivated and managed vegetation / agriculture' = 40,
  'Urban / built up' = 50,
  'Bare / sparse vegetation' = 60,
  'Snow and ice' = 70,
  'Permanent water bodies' = 80,
  'Herbaceous wetland' = 90,
  'Moss and lichen' = 100,
  'Closed forest, evergreen needle leaf' = 111,
  'Closed forest, evergreen broad leaf' = 112,
  'Closed forest, deciduous needle leaf' = 113,
  'Closed forest, deciduous broad leaf' = 114,
  'Closed forest, mixed' = 115,
  'Closed forest, not matching any of the other definitions' = 116,
  'Open forest, evergreen needle leaf' = 121,
  'Open forest, evergreen broad leaf' = 122,
  'Open forest, deciduous needle leaf' = 123,
  'Open forest, deciduous broad leaf' = 124,
  'Open forest, mixed' = 125,
  'Open forest, not matching any of the other definitions' = 126,
  'Oceans, seas' = 200
)

# Get coverage radar
get_coverage_radar <- function(radar_file) {
  sum(is.na(valuesRadar) + 1)/(511*511)
}

# Detect s1 data
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

#' Download thumbnails
download_thumbnails <- function(thumbnail_file) {
  in_thumb_file <- paste0(tempfile(), ".tif")
  out_thumb_file <- paste0(tempfile(), ".tif")
  drive_download(
    file = thumbnail_file,
    path = in_thumb_file,
    overwrite = TRUE
  )
  system(sprintf("gdalwarp  %s %s -of COG -co BLOCKSIZE=256", in_thumb_file, out_thumb_file))
  out_thumb_file
}


#' Download and order and save cloud models in the temp/ folder
files_target_final_creator <- function(files_target_files, s2_id, point) {
  # s2cloudness_prob
  s2cloudness_file <- sprintf("%s.tif", tempfile())
  from_TIFF_to_COG(files_target_files[["s2cloudness_prob"]], s2cloudness_file)

  # sen2cor_real
  sen2cor_file <- sprintf("%s.tif", tempfile())
  from_TIFF_to_COG(files_target_files[["sen2cor_real"]], sen2cor_file)

  # manual
  manual_file <- sprintf("%s.tif", tempfile())
  from_TIFF_to_COG(files_target_files[["manual"]], manual_file)

  # IPLcloud
  inIPLcloud <- sprintf("%s.tif", tempfile())
  outIPLcloud <- sprintf("%s.tif", tempfile())
  qa_band <- raster(files_target_files[["IPL_cloudmask_reclass"]])
  qa_band[1,1] <- 1
  qa_band[1,3] <- 3
  qa_band[qa_band==3] <- 1
  qa_band[qa_band==1] <- 2
  writeRaster(qa_band, inIPLcloud, overwrite = TRUE)
  from_TIFF_to_COG(inIPLcloud, outIPLcloud)

  # QAbands
  inQAbands <- sprintf("%s.tif", tempfile())
  outQAbands <- sprintf("%s.tif", tempfile())
  rK <- 1000
  sen2_q60 <- ee$Image(sprintf("COPERNICUS/S2/%s", s2_id))$select("QA60")
  cloud_mask <- sen2_q60$remap(c(0L, 1024L, 2048L), c(0L, 1L, 2L))
  new_bound <- (st_bbox(read_stars(outIPLcloud)) + c(-rK, -rK, rK, rK)) %>% st_as_sfc()
  rr1 <- ee_as_raster(
    image = cloud_mask,
    region = new_bound %>% sf_as_ee()
  )
  tempfiletif <- sprintf("%s.tif", tempfile())
  writeRaster(rr1, tempfiletif)
  system(
    sprintf(
      "gdalwarp %s %s -overwrite -te %s %s %s %s -tr 10 10",
      tempfiletif,
      inQAbands,
      as.numeric(st_bbox(read_stars(outIPLcloud)))[1],
      as.numeric(st_bbox(read_stars(outIPLcloud)))[2],
      as.numeric(st_bbox(read_stars(outIPLcloud)))[3],
      as.numeric(st_bbox(read_stars(outIPLcloud)))[4]
    )
  )
  from_TIFF_to_COG(inQAbands, outQAbands)

  # FMASK (csaybar folder)
  inFMASK <- sprintf("%s.tif", tempfile())
  outFMASK <- sprintf("%s.tif", tempfile())

  files_points_general <- googledrive::drive_ls(
    path = as_id("1DV7jPg5jyxc58GdbK18LueQ1c_ksHxem"),
    q = sprintf("name contains '%s'", point)
  )
  files_s2_files <- googledrive::drive_ls(
    path = as_id(files_points_general$id)
  )
  drive_download(
    file = files_s2_files[grepl(s2_id, files_s2_files$name), ],
    path = inFMASK
  )
  from_TIFF_to_COG(inFMASK, outFMASK)

  # SIAM
  inSIAM <- sprintf("%s.tif", tempfile())
  outSIAM <- sprintf("%s.tif", tempfile())

  files_points_general <- googledrive::drive_ls(
    path = as_id("1DXDzgX91PkGjYXg5HZfE1QEeDaVjm6kU"),
    q = sprintf("name contains '%s'", point)
  )
  files_s2_files <- googledrive::drive_ls(
    path = as_id(files_points_general$id)
  )
  drive_download(
    file = files_s2_files[grepl(s2_id, files_s2_files$name), ],
    path = inSIAM
  )
  from_TIFF_to_COG(inSIAM, outSIAM)

  list(
    manual = manual_file,
    s2cloudness = s2cloudness_file,
    sen2cor = sen2cor_file,
    iplcloud = outIPLcloud,
    qa60 = outQAbands,
    fmask = outFMASK,
    siam = outSIAM
  )
}

#' From GeoTIFF to COG
from_TIFF_to_COG <- function(infile, outfile) {
  system(sprintf("gdalwarp  %s %s -of COG -overwrite -co BLOCKSIZE=256", infile, outfile))
  outfile
}

#' From GeoTIFF to COG
create_final_target <- function(x, output) {
  clear <- x == 0
  thick_cloud <- x == 1
  thin_cloud <- x == 2
  shadow_cloud <- x == 3
  stk_stars <- stack(clear, thick_cloud, thin_cloud, shadow_cloud) %>%
    st_as_stars() %>%
    '[['("layer.1") %>%
    np$save(output, .)
  output
}


#' Create .npy files from input model (20 GeoTIFF files)
#' @noRd
create_npy_file <- function(files_input_image_ls, output) {
  # Names order
  names_in_order <- c(
    "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11",
    "B12", "VH", "VV", "angle", "CDI", "cloudshadow_direction", "elevation",
    "landuse"
  )
  files_input_image_ls[names_in_order] %>%
    read_stars() %>%
    merge() %>%
    '[['("X") %>%
    np$save(output, .)
  output
}

#' Download input from Roy Folder
download_input <- function(files_input_image) {
  list_manual <- list()
  for (index in seq_len(nrow(files_input_image))) {
    file_tif <- paste0(tempfile(), ".tif")
    drive_download(
      file = files_input_image[index,],
      path = file_tif,
      overwrite = TRUE
    )
    fn_name <- gsub("\\.tif$", "", files_input_image[index,]$name)
    list_manual[[fn_name]] <- file_tif
  }
  list_manual
}

#' List images in a Point Folder
list_point_folder <- function(point_name, gd_id = "1BeVp0i-dGSuBqCQgdGZVDj4qzX1ms7L6") {

  # refresh
  httr::set_config(httr::config( ssl_verifypeer = 0L))
  httr::set_config(httr::config(http_version = 0))

  # cloudsen12 level find the point
  files_points_general <- googledrive::drive_ls(
    path = as_id(gd_id),
    q = sprintf("name contains '%s'", point_name)
  )

  # POINT level
  googledrive::drive_ls(
    path = as_id(files_points_general$id)
  )
}

#' List images in a S2 Image Folder
list_images_in_folder <- function(point_name) {
  lpf <- list_point_folder(point_name)
  lpf[grepl("^[0-9]", lpf$name),]
}

#' List a specific image
list_image_folder <- function(image_id) {
  googledrive::drive_ls(
    path = as_id(image_id$id)
  )
}

# List a specific image input folder
list_input_image_folder <- list_image_folder

# List a specific image target folder
list_target_image_folder <- list_image_folder

# List a specific image thumbnails folder
list_thumbnails_image_folder <- list_image_folder

#' Download cloud model results :)
download_labels <- function(files_target_image) {
  list_manual <- list()
  for (index in seq_len(nrow(files_target_image))) {
    file_tif <- paste0(tempfile(), ".tif")
    drive_download(
      file = files_target_image[index,],
      path = file_tif,
      overwrite = TRUE
    )
    fn_name <- gsub("\\.tif$", "", files_target_image[index,]$name)
    list_manual[[fn_name]] <- file_tif
  }

  if (is.null(list_manual[["IPL_cloudmask_reclass"]])) {
    baseRaster <- raster(list_manual[["manual"]])
    baseRaster[] <- NA
    ipl_temp <- paste0(tempfile(), ".tif")
    writeRaster(baseRaster, ipl_temp)
    list_manual[["IPL_cloudmask_reclass"]] <- ipl_temp
  }
  list_manual
}

#' DataBase Migration
db_migration <- function(point, output) {

  # Point level (ROY FOLDER)
  files_points <- list_point_folder(point_name = tolower(point))
  files_points_only_img_folder <- files_points[grepl("^[0-9]", files_points$name),]

  # List each image
  for (index in 1:5) {
    # Image level (Roy Folder)
    s2_id <- files_points_only_img_folder[index,]$name

    # List files inside the IMAGE Folder
    files_image_folder <- list_image_folder(files_points_only_img_folder[index,])

    # Target level (Roy Folder)
    files_target_folder <- list_target_image_folder(files_image_folder[files_image_folder$name == "target",])

    # Download target image (Roy Folder)
    files_target_files <- download_labels(files_target_folder)
    manual_file <- files_target_files$manual
    target_npy_temp <- create_final_target(raster(manual_file), paste0(tempfile(), ".npy"))
    target_model_folder <- files_target_final_creator(files_target_files, s2_id, point)

    # Input level (Roy Folder)
    files_input_folder <- list_target_image_folder(files_image_folder[files_image_folder$name == "input",])
    files_input_folder_no_b10_thershold <- files_input_folder %>%
      filter(name != "B10_threshold.tif")
    files_input_files <- download_input(files_input_folder_no_b10_thershold)
    input_npy_temp <- create_npy_file(files_input_files, paste0(tempfile(), ".npy")) # Creat npy array (20x511x511)

    # Thumbnails level (Roy Folder)
    thumbnails_image_list <- list_thumbnails_image_folder(files_image_folder[files_image_folder$name == "thumbnails",])
    thumbnail_file <- thumbnails_image_list[thumbnails_image_list$name == "thumbnail.tif",]
    thumbnail_tif_temp <- download_thumbnails(thumbnail_file)

    # Create a directory
    dir.create(
      path = sprintf("%s/%s/%s/models",output, point, s2_id),
      recursive = TRUE,
      showWarnings = FALSE
    )

    # save thumbnail
    file.copy(
      from = thumbnail_tif_temp,
      to = sprintf("%s/%s/%s/thumbnail.tif",output, point, s2_id),
      overwrite = TRUE
    )

    # save manual npy
    file.copy(
      from = target_npy_temp,
      to = sprintf("%s/%s/%s/manual.npy",output, point, s2_id),
      overwrite = TRUE
    )

    # save input npy
    file.copy(
      from = input_npy_temp,
      to = sprintf("%s/%s/%s/input.npy",output, point, s2_id),
      overwrite = TRUE
    )

    # models
    for (index in seq_along(target_model_folder)) {
      # save input npy
      file.copy(
        from = target_model_folder[[index]],
        to = sprintf("%s/%s/%s/models/%s.tif",output, point, s2_id, names(target_model_folder[index])),
        overwrite = TRUE
      )
    }

    # Point Level (Cesar Folder)
    # files_points_general <- googledrive::drive_ls(
    #   path = as_id("1e4C8pWF6GiLOchzUg5urjptSS8vWVxRH"),
    #   q = sprintf("name contains '%s'", point)
    # )
    #
    # # Create a new folder with the image ID
    # drive_auth("s1078735@stud.sbg.ac.at")
    # new_folder <- drive_mkdir(
    #   path = files_points_general,
    #   name = s2_id,
    #   overwrite = TRUE
    # )
    #
    # # Upload thumbnail to CLOUDSEN12
    # drive_upload(
    #   media = thumbnail_tif_temp,
    #   path = as_id(new_folder),
    #   name = "thumbnail.tif",
    #   overwrite = TRUE
    # )
    #
    # # Upload manual target to CLOUDSEN12
    # drive_upload(
    #   media = target_npy_temp,
    #   path = as_id(new_folder),
    #   name = "manual.npy",
    #   overwrite = TRUE
    # )
    #
    # # Upload input to CLOUDSEN12
    # drive_upload(
    #   media = input_npy_temp,
    #   path = as_id(new_folder),
    #   name = "input.npy",
    #   overwrite = TRUE
    # )
    #
    # # Create model folder
    # models_folder <- drive_mkdir(
    #   path = new_folder,
    #   name = "models",
    #   overwrite = TRUE
    # )
    #
    # # Upload cloud models
    # for (index in seq_along(target_model_folder)) {
    #   models_link <- drive_upload(
    #     media = target_model_folder[[index]],
    #     path = as_id(models_folder$id),
    #     name = paste0(names(target_model_folder[index]), ".tif"),
    #     overwrite = TRUE
    #   )
    #   drive_share_anyone(models_link, verbose = FALSE)
    # }
  }
}

# STAC ITEM creator
stac_feature_creator <- function(point,  potencial_points, output = output) {
  # Get the number i.e Point_1000 -> 1000
  pt_number <- as.numeric(strsplit(point, "_")[[1]][2])

  # Read XLS  (columns : point, labeler, 	type, 	exist, 	difficulty, sen2_id) (ROY FOLDER)
  # https://docs.google.com/spreadsheets/d/1LpW9JY2BdhlQvAObD1BCzoBiliNRnU3fCRWMQJFRvoM
  start_read <- (5*(pt_number - 1) + 2)
  sheet_range <- sprintf("A%s:F%s", start_read, (start_read + 4))
  header <-  c("point", "labeler", "type", "exist", "difficulty", "sen2_id")
  sheet_xls <- googlesheets4::read_sheet(
    ss = "1LpW9JY2BdhlQvAObD1BCzoBiliNRnU3fCRWMQJFRvoM",
    range = sheet_range,
    col_names = header
  )

  # Point level (ROY FOLDER)
  # https://drive.google.com/drive/u/1/folders/1aTPIZ974zvtti6a02eiMyZaIf_Rp8QEc
  files_points <- list_point_folder(point)
  files_points_only_img_folder <- files_points[grepl("^[0-9]", files_points$name),]

  for (index in 1:5) {

    ### GENEARAL PROPERTIES ---------------------

    # sentinel2_product_id_gee
    sentinel2_product_id_gee <- files_points_only_img_folder$name[index]

    # sentinel2_product_id
    s2_id <- ee$Image(sprintf("COPERNICUS/S2/%s", sentinel2_product_id_gee))
    sentinel2_product_id <- s2_id$get("PRODUCT_ID")$getInfo()

    # sentinel2_date
    sentinel2_date <- ee_get_date_img(s2_id)[["time_start"]]
    sentinel2_date_f <- paste0(
      paste(strsplit(sentinel2_date %>% as.character(), " ")[[1]], collapse = "T"),
      "Z"
    )

    # sentinel1_product_id
    metadata_name <- paste0("metadata_", strsplit(point, "_")[[1]][2], ".json")
    json_file <- googledrive::drive_ls(
      path = as_id("1fBGAjZkjPEpPr0p7c-LtJmfbLq3s87RK"),
      q = sprintf("name contains '%s'", metadata_name)
    ) %>% googledrive::drive_download(path = tempfile())
    jsonfile_r <- jsonlite::read_json(json_file$local_path)

    ## geom
    st_point <- st_sfc(geometry = st_point(c(jsonfile_r$x, jsonfile_r$y)), crs = 4326)
    crs_kernel <- ee$Image(s2_id)$select(0)$projection()$getInfo()$crs
    point_utm <- st_transform(st_point, crs_kernel)
    ee_point <- ee$Geometry$Point(point_utm[[1]], proj = crs_kernel)

    ## sentinel1_product_id
    sentinel1_product_id <- ee_get_s1(point = ee_point, s2_date = sentinel2_date)
    sentinel1_date <- ee_get_date_img(ee$Image(sentinel1_product_id))[["time_start"]]
    sentinel1_date_f <- paste0(
      paste(strsplit(sentinel1_date %>% as.character(), " ")[[1]], collapse = "T"),
      "Z"
    )

    # label_type
    label_type <- potencial_points[pt_number,]$label

    # difficulty
    difficulty <- sheet_xls %>%
      dplyr::filter(sen2_id==sentinel2_product_id_gee) %>%
      '[['('difficulty')

    # annotator_name
    annotator_name <- sheet_xls %>%
      dplyr::filter(sen2_id==sentinel2_product_id_gee) %>%
      '[['('labeler')

    # grd_post_processing_software_name
    grd_post_processing_software_name <- ee$Image(sentinel1_product_id)$
      get("GRD_Post_Processing_software_name")$
      getInfo()

    # grd_post_processing_software_version
    grd_post_processing_software_version <- ee$Image(sentinel1_product_id)$
      get("GRD_Post_Processing_software_version")$
      getInfo()

    # slc_processing_facility_name
    slc_processing_facility_name <- ee$Image(sentinel1_product_id)$
      get("SLC_Processing_facility_name")$
      getInfo()

    # SLC_Processing_software_version
    slc_processing_software_version <- ee$Image(sentinel1_product_id)$
      get("SLC_Processing_software_version")$
      getInfo()

    # fmask_version
    fmask_version <- "4.3.0"

    # sencloudness_version
    sencloudness_version <- "1.5.0"

    # sen2cor_version
    ee_s2sr <- ee$Image(sprintf("COPERNICUS/S2_SR/%s", sentinel2_product_id_gee))
    sen2cor_dtidentf <- ee_s2sr$get("DATATAKE_IDENTIFIER")$getInfo()
    sen2cor_version <- strsplit(sen2cor_dtidentf, "_")[[1]][4]

    # iplcloud_version
    iplcloud_version <- "0.1"

    # iplcloud_nimg
    iplcloud_nimg <- 3

    # iplcloud_method
    iplcloud_method <- "persistence"

    # cloud_thickness
    cloud_thickness <- jsonfile_r[[sentinel2_product_id_gee]][["cloud_thickness"]]

    # radar_coverage
    input_data <- np$load(sprintf("%s/points/%s/%s/input.npy", dirname(output), point, sentinel2_product_id_gee))
    radar_coverage <- sum(!is.na(input_data[,,14]))/(511*511) * 100

    # cloud_height
    cloud_height <- jsonfile_r[[sentinel2_product_id_gee]][["cloud_height"]]

    # cloud_type
    cloud_type <- jsonfile_r[[sentinel2_product_id_gee]][["cloud_type"]]

    ### VIEW extention ---------------------

    # view_off_nadir
    view_off_nadir <- 0

    # view_sun_azimuth
    view_sun_azimuth <- s2_id$get("MEAN_SOLAR_AZIMUTH_ANGLE")$getInfo()

    # view_sun_elevation
    view_sun_elevation <- 90 - s2_id$get("MEAN_SOLAR_ZENITH_ANGLE")$getInfo()

    ### PROJ extention ---------------------
    raster_base_file <- sprintf(
      "%s/points/%s/%s/models/fmask.tif",
      dirname(output), point, sentinel2_product_id_gee
    )

    sf_tile_file <- tempfile()
    sf_tile <- suppressWarnings(
      drive_download(
        files_points[grepl("\\.gpkg$", files_points$name), ],
        sf_tile_file, overwrite = TRUE,
      )$local_path %>% st_read()
    )

    # proj_shape
    proj_shape <- c(511, 511)

    # proj_epsg
    proj_epsg <- st_crs(sf_tile)$epsg

    # proj_geometry
    proj_geometry <- geojsonio::geojson_list(sf_tile$geom)
    attributes(proj_geometry) <- list(names = c("type", "coordinates"))

    # proj_centroid
    df_xy <- st_coordinates(st_transform(st_centroid(sf_tile), 4326)[["geom"]]) %>%
      as.data.frame()
    proj_centroid <- list(
      "lat" = df_xy[["Y"]],
      "lon" = df_xy[["X"]]
    )

    # proj_transform
    proj_transform <- ee_utils_py_to_r(rio$open(raster_base_file)$get_transform())


    ### CLOUD COVERAGE PROPERTIES ---------------------
    manual_base_file <- sprintf(
      "%s/points/%s/%s/manual.npy",
      dirname(output), point, sentinel2_product_id_gee
    )

    if (file.exists(manual_base_file) & (label_type == "high_quality")) {
      manual_labelling <- np$load(manual_base_file)

      # clear_mean_coverage
      clear_mean_coverage <- sum(manual_labelling[,,1])/(511*511)

      # thickcloud_mean_coverage
      thickcloud_mean_coverage <- sum(manual_labelling[,,2])/(511*511)

      # thincloud_cloud_coverage
      thincloud_cloud_coverage <- sum(manual_labelling[,,3])/(511*511)

      # cloudshadow_mean_coverage
      cloudshadow_mean_coverage <- sum(manual_labelling[,,4])/(511*511)

    } else {
      # clear_mean_coverage
      clear_mean_coverage <- NA

      # thickcloud_mean_coverage
      thickcloud_mean_coverage <- NA

      # thincloud_cloud_coverage
      thincloud_cloud_coverage <- NA

      # cloudshadow_mean_coverage
      cloudshadow_mean_coverage <- NA
    }

    ### MEAN VALUES PROPERTIES ---------------------
    input <- np$load(
      sprintf(
        "%s/points/%s/%s/input.npy",
        dirname(output), point, sentinel2_product_id_gee
      )
    )

    # B1
    b1_mean <- mean(input[,,1], na.rm = TRUE)

    # B2
    b2_mean <- mean(input[,,2], na.rm = TRUE)

    # B3
    b3_mean <- mean(input[,,3], na.rm = TRUE)

    # B4
    b4_mean <- mean(input[,,4], na.rm = TRUE)

    # B5
    b5_mean <- mean(input[,,5], na.rm = TRUE)

    # B6
    b6_mean <- mean(input[,,6], na.rm = TRUE)

    # B7
    b7_mean <- mean(input[,,7], na.rm = TRUE)

    # B8
    b8_mean <- mean(input[,,8], na.rm = TRUE)

    # B8A
    b8a_mean <- mean(input[,,9], na.rm = TRUE)

    # B9
    b9_mean <- mean(input[,,10], na.rm = TRUE)

    # B10
    b10_mean <- mean(input[,,11], na.rm = TRUE)

    # B11
    b11_mean <- mean(input[,,12], na.rm = TRUE)

    # B12
    b12_mean <- mean(input[,,13], na.rm = TRUE)

    # VH
    vh_mean <- mean(input[,,14], na.rm = TRUE)

    # VV
    vv_mean <- mean(input[,,15], na.rm = TRUE)

    # angle
    angle_mean <- mean(input[,,16], na.rm = TRUE)

    # CDI
    cdi_mean <- mean(input[,,17], na.rm = TRUE)

    # cloudshadow_direction
    cloudshadow_direction <- mean(input[,,18], na.rm = TRUE)

    # elevation
    elevation_mean <- mean(input[,,19], na.rm = TRUE)

    # landuse
    # land_use
    landuse_mode <- Modes(as.numeric(input[,,20]))
    land_use <- names(landuse_types_list[landuse_types_list == landuse_mode])

    final_json_properties <- list(
      point_id = point,
      sentinel2_product_id = sentinel2_product_id,
      sentinel2_product_id_gee = sentinel2_product_id_gee,
      sentinel2_date = sentinel2_date_f,
      datetime = sentinel2_date_f,
      sentinel1_product_id = basename(sentinel1_product_id),
      sentinel1_date = sentinel1_date_f,
      label_type = label_type,
      difficulty = difficulty,
      annotator_name = annotator_name,
      fmask_version = fmask_version,
      sen2cor_version = sen2cor_version,
      sencloudness_version = sencloudness_version,
      iplcloud_version = iplcloud_version,
      iplcloud_nimg = iplcloud_nimg,
      iplcloud_method = iplcloud_method,
      slc_processing_facility_name = slc_processing_facility_name,
      slc_processing_software_version = slc_processing_software_version,
      grd_post_processing_software_name = grd_post_processing_software_name,
      grd_post_processing_software_version = grd_post_processing_software_version,
      radar_coverage = radar_coverage,
      cloud_height = cloud_height,
      cloud_type = cloud_type,
      land_use = land_use,
      'view:off_nadir' = view_off_nadir,
      'view:sun_azimuth' = view_sun_azimuth,
      'view:sun_elevation' = view_sun_elevation,
      'proj:epsg' = proj_epsg,
      'proj:geometry' = proj_geometry,
      'proj:shape' = proj_shape,
      'proj:centroid' = proj_centroid,
      'proj:transform' = proj_transform,
      b1_mean = b1_mean,
      b2_mean = b2_mean,
      b3_mean = b3_mean,
      b4_mean = b4_mean,
      b5_mean = b5_mean,
      b6_mean = b6_mean,
      b7_mean = b7_mean,
      b8_mean = b8_mean,
      b8a_mean = b8a_mean,
      b9_mean = b9_mean,
      b10_mean = b10_mean,
      b11_mean = b11_mean,
      b12_mean = b12_mean,
      vh_mean = vh_mean,
      vv_mean = vv_mean,
      angle_mean = angle_mean,
      cdi_mean = cdi_mean,
      cloudshadow_direction = cloudshadow_direction,
      elevation_mean = elevation_mean,
      landuse_mode = landuse_mode,
      clear_mean_coverage = clear_mean_coverage,
      thickcloud_mean_coverage = thickcloud_mean_coverage,
      thincloud_cloud_coverage = thincloud_cloud_coverage,
      thincloud_cloud_coverage = cloudshadow_mean_coverage
    )

    # ID list
    id_f <- paste(point, sentinel2_product_id_gee,sep = "_")

    # geometry WGS84
    geometry_f <- st_transform(sf_tile$geom, 4326) %>%
      geojsonio::geojson_list()
    attributes(geometry_f) <- list(names = c("type", "coordinates"))

    # bbox
    bbox_f <- st_transform(sf_tile$geom, 4326) %>% st_bbox()
    attributes(bbox_f) <- NULL


    # Remove ---------
    local_folder_files <- sprintf("%s/points/%s/%s/", dirname(output), point, sentinel2_product_id_gee)

    href_thumbnail <- paste0(local_folder_files, "thumbnail.tif")
    href_fmask <- paste0(local_folder_files, "/models/fmask.tif")
    href_iplcloud <- paste0(local_folder_files, "/models/iplcloud.tif")
    href_manual <- paste0(local_folder_files, "/models/manual.tif")
    href_qa60 <- paste0(local_folder_files, "/models/qa60.tif")
    href_s2cloudness <- paste0(local_folder_files, "/models/s2cloudness.tif")
    href_sen2cor <- paste0(local_folder_files, "/models/sen2cor.tif")
    href_siam <- paste0(local_folder_files, "/models/siam.tif")

    final_json <- list(
      id = id_f,
      type = "Feature",
      collection = "cloudsen12",
      links = list(),
      geometry = geometry_f,
      properties = final_json_properties,
      assets = list(
        thumbnail = list(
          type = "image/tiff; application=geotiff; profile=cloud-optimized",
          href = href_thumbnail,
          title = "Thumbnail image - COG"
        ),
        fmask = list(
          type = "image/tiff; application=geotiff; profile=cloud-optimized",
          href = href_fmask,
          title = "Fmask label - COG"
        ),
        iplcloud = list(
          type = "image/tiff; application=geotiff; profile=cloud-optimized",
          href = href_iplcloud,
          title = "IPLcloud label - COG"
        ),
        manual = list(
          type = "image/tiff; application=geotiff; profile=cloud-optimized",
          href = href_manual,
          title = "manual label - COG"
        ),
        qa60 = list(
          type = "image/tiff; application=geotiff; profile=cloud-optimized",
          href = href_qa60,
          title = "qa60 label - COG"
        ),
        s2cloudness = list(
          type = "image/tiff; application=geotiff; profile=cloud-optimized",
          href = href_s2cloudness,
          title = "s2cloudness label - COG"
        ),
        sen2cor = list(
          type = "image/tiff; application=geotiff; profile=cloud-optimized",
          href = href_sen2cor,
          title = "sen2cor label - COG"
        ),
        siam = list(
          type = "image/tiff; application=geotiff; profile=cloud-optimized",
          href = href_siam,
          title = "siam label - COG"
        )
      ),
      bbox = bbox_f,
      stac_extensions = list("proj", "view"),
      stac_version = "1.0.0-beta.2"
    )

    # Create a directory
    folder_f <- sprintf("%s/%s/",output, point)
    dir.create(
      path = folder_f,
      recursive = TRUE,
      showWarnings = FALSE
    )

    # Write final JSON :)
    jsonlite::write_json(
      x = final_json,
      path = sprintf("%s/%s.json", folder_f, sentinel2_product_id_gee),
      pretty = TRUE,
      auto_unbox = TRUE
    )
  }
}

#' Upgrade asset in Item STAC
upgrade_link_drive <- function(point, output) {
  # sentinel1_product_id
  metadata_name <- paste0("metadata_", strsplit(point, "_")[[1]][2], ".json")
  json_file <- googledrive::drive_ls(
    path = as_id("1fBGAjZkjPEpPr0p7c-LtJmfbLq3s87RK"),
    q = sprintf("name contains '%s'", metadata_name)
  ) %>% googledrive::drive_download(path = tempfile())
  jsonfile_r <- jsonlite::read_json(json_file$local_path)
  s2_ids <- names(jsonfile_r[grepl("^[0-9]", names(jsonfile_r))])

  # Folder to find Public links (CSAYBAR FOLDER)
  files_points_general <- googledrive::drive_ls(
    path = as_id("1tcFgbP3SLovBy3UNs7OPMWNOO5PuiCuP"),
    q = sprintf("name contains '%s'", point)
  )
  files_specific_folder <- googledrive::drive_ls(
    path = files_points_general
  )


  index <- 3
  for (index in 1:5) {
    jsfile <- jsonlite::read_json(
      path = sprintf("%s/jsons/%s/%s.json",output, point, s2_ids[index])
    )

    # assets
    image_in_drive <- files_specific_folder[files_specific_folder$name %in% s2_ids[index],]
    image_files <- drive_ls(image_in_drive, recursive = TRUE)

    # Remove ---------
    href_thumbnail_img <- image_files[image_files$name %in% "thumbnail.tif",]
    href_thumbnail <- sprintf("https://drive.google.com/uc?id=%s&export=download", href_thumbnail_img$id)

    href_fmask_img <- image_files[image_files$name %in% "fmask.tif",]
    href_fmask <- sprintf("https://drive.google.com/uc?id=%s&export=download", href_fmask_img$id)

    href_iplcloud_img <- image_files[image_files$name %in% "iplcloud.tif",]
    href_iplcloud <- sprintf("https://drive.google.com/uc?id=%s&export=download", href_iplcloud_img$id)

    href_manual_img <- image_files[image_files$name %in% "manual.tif",]
    href_manual <- sprintf("https://drive.google.com/uc?id=%s&export=download", href_manual_img$id)

    href_qa60_img <-  image_files[image_files$name %in% "qa60.tif",]
    href_qa60 <- sprintf("https://drive.google.com/uc?id=%s&export=download", href_qa60_img$id)

    href_s2cloudness_img <- image_files[image_files$name %in% "s2cloudness.tif",]
    href_s2cloudness <- sprintf("https://drive.google.com/uc?id=%s&export=download", href_s2cloudness_img$id)

    href_sen2cor_img <- image_files[image_files$name %in% "sen2cor.tif",]
    href_sen2cor <- sprintf("https://drive.google.com/uc?id=%s&export=download", href_sen2cor_img$id)

    href_siam_img <-  image_files[image_files$name %in% "siam.tif",]
    href_siam <- sprintf("https://drive.google.com/uc?id=%s&export=download", href_siam_img$id)

    # upgrade json
    jsfile$assets$thumbnail$href <- href_thumbnail
    jsfile$assets$fmask$href <- href_fmask
    jsfile$assets$iplcloud$href <- href_iplcloud
    jsfile$assets$manual$href <- href_manual
    jsfile$assets$qa60$href <- href_qa60
    jsfile$assets$s2cloudness$href <- href_s2cloudness
    jsfile$assets$sen2cor$href <- href_sen2cor
    jsfile$assets$siam$href <- href_siam

    # Write final JSON :)
    jsonlite::write_json(
      x = jsfile,
      path = sprintf("%s/jsons/%s/%s.json",output, point, s2_ids[index]),
      pretty = TRUE,
      auto_unbox = TRUE
    )
  }
}


# FeatureCollection Creator
stac_featurecollection_creator <- function(json_list, outputfile) {
  fc_f <- list(
    type = "FeatureCollection",
    features = lapply(json_list, jsonlite::read_json)
  )
  # Write final JSON :)
  jsonlite::write_json(
    x = fc_f,
    path = outputfile,
    pretty = TRUE,
    auto_unbox = TRUE
  )
}


# DB migration batch
db_migration_batch <- function(points, output) {
  # Create folder
  if (!dir.exists(output)) {
    dir.create(path = paste0(output, "/points"), recursive = TRUE, showWarnings = TRUE)
  }
  for (index in points) {
    message(sprintf("Working in the Point_%04d :D", index))
    point <- sprintf("Point_%04d", index)
    try(db_migration(point = point, output = sprintf("%s/%s", output, "points")))
  }
}

# STAC features creator
stac_feature_creator_batch <- function(points, potencial_points, output) {
  output <- sprintf("%s/%s", output, "jsons")
  # Create folder
  if (!dir.exists(output)) {
    dir.create(path = output, recursive = TRUE, showWarnings = TRUE)
  }

  for (index in points) {
    message(sprintf("Working with STAC :D --- Point_%04d", index))
    point <- sprintf("Point_%04d", index)
    try(
      stac_feature_creator(
        point = point,
        potencial_points = potencial_points,
        output = output
      )
    )
  }
}
