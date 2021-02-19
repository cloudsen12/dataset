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
      ) %>% ee$Image$unmask(-99)
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
      lab_lab_user(path = dirname(metadata_main), point_name = point_name)
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
  s1_grd <- ee$Image(s1_id)$unmask(-99)
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
  s2_img432 <- ee_as_raster(
    image = ee$Image(s2_id)$select(c("B4","B3","B2")),
    scale = 500,
    quiet = TRUE
  )
  roi <- extent(final_stack[[1]]) %>%
    st_bbox() %>%
    st_as_sfc()
  st_crs(roi) <- crs_kernel
  png(sprintf("%s/thumbnails/thumbnail.png", output_final_folder), 1000, 1000)
  max_value <- max(maxValue(s2_img432))
  plotRGB(s2_img432/max_value, r = 3, g = 2, b = 1, stretch = "lin")
  plot(roi, add=TRUE, border = "red", lwd = 3)
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
upload_results <- function(folder_id, pnt_name) {

  # files Google Drive ID
  files_points <- googledrive::drive_ls(
    path = as_id(folder_id)
  )

  # List files (Get SENTINEL_2 ID)
  directories <- list.files()
  directories_sen2ids <- directories[!grepl("iris|viz", directories)]
  directories_sen2ids <- directories_sen2ids[dir.exists(directories_sen2ids)]

  # Get the point number
  point_name <- gsub("\\.gpkg$", "", directories[grepl("gpkg", directories)])

  if (point_name  != pnt_name) {
    stop("The Google Drive folder ID doesn't match with the point name.")
  }

  # Reference raster (base_raster)
  reference_raster <- raster(paste0(directories_sen2ids[1],"/input/B1.tif"))


  for (directory in directories_sen2ids) {

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

#' Generate R script
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
      "",
      "source(\"https://gist.githubusercontent.com/csaybar/8a4487f1fd1c488be3dce0b60f7f0ce8/raw/5353cf8f5cc0bb5a76375b4895bd3ad4dd71b986/cloudsen12_functions.R\")",
      "",
      "# Generate svg for each image",
      "generate_preview()",
      "",
      "# Upload your results to Google Drive",
      "upload_results()"
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
