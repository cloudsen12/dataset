#' CloudSEN12 thumbnail creator
#'
#' Create several thumbnail to find the clouds closest to the theoretical probabilities
#' All Sentinel2 images have a Sentinel1 pair with no more than 2.5 days of delay.
#'
select_dataset_thumbnail_creator <- function(cloudsen2_row,
                                             n_images = 50,
                                             kernel_size = c(255, 255),
                                             data_range = c("2019-01-01", "2020-07-31"),
                                             output = "results/") {
  # 1. Create output directory
  dir.create(output, showWarnings = FALSE)

  # 2. Create a point which represent the center of the chip (from local to ee)
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


metadata_dataset_creator <- function(cloudsen2_row,
                                     output = "results/") {
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

dataset_creator_chips <- function(jsonfile,
                                  kernel_size = c(255, 255),
                                  output_final = "final_results/") {
  # 1. Read JSON file
  jsonfile_r <- jsonlite::read_json(jsonfile)

  # 2. Create a point which represent the center of the chip
  st_point <- st_sfc(geometry = st_point(c(jsonfile_r$x, jsonfile_r$y)), crs = 4326)
  point <- ee$Geometry$Point(jsonfile_r$x, jsonfile_r$y)

  # 3. Identify all the S2 images
  s2_ids <- sprintf("COPERNICUS/S2/%s", names(jsonfile_r)[1:5])

  # s2_id <- s2_ids[2]
  # 3. Download each image at each point
  for (s2_id in s2_ids) {
    message(sprintf("Downloading: %s", s2_id))
    # 3.1 S2 ID and dates
    s2_img <- ee$Image(s2_id)
    # Map$centerObject(point)
    # map02 <- Map$addLayer(s2_img, list(min=0, max=10000, bands = c("B4","B3","B2"))) +
    # Map$addLayer(point)
    s2_date <- ee_get_date_img(s2_img)[["time_start"]]

    # 3.2 S1 ID and dates
    s1_id <- ee_get_s1(point = point, s2_date = s2_date)
    s1_img <- ee$Image(s1_id)
    # Map$addLayer(s1_img)

    # 3.3 Create a S2 Image with cloud mask information
    s2_fullinfo <- ee_merge_s2_full(s2_id, s1_id, s2_date)
    crs_kernel <- s2_fullinfo$select(0)$projection()$getInfo()$crs
    point_utm <- st_transform(st_point, crs_kernel)
    ee_point <- ee$Geometry$Point(point_utm[[1]], proj = crs_kernel)

    # 3.4 Create a 511x511 tile
    band_names <- c(s2_fullinfo$bandNames()$getInfo(), "x", "y")
    s2_img_array <- s2_fullinfo$addBands(s1_img) %>%
      ee$Image$addBands(ee$Image$pixelCoordinates(projection = crs_kernel)) %>%
      ee$Image$neighborhoodToArray(
        kernel = ee$Kernel$rectangle(kernel_size[1], kernel_size[2], "pixels")
      ) %>%
      ee$Image$sampleRegions(ee$FeatureCollection(point),
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

    ### Prepare data for iris ------------------------
    output_final_d <- sprintf("%s/images", output_final)
    output_final_folder <- sprintf("%s/images/%s", output_final, basename(s2_id))
    metadata_final <- sprintf("%s/cloud-segmentation.json", output_final)

    metadata_spec <- sprintf("%s/images/%s/metadata.json", output_final, basename(s2_id))
    inputdata_spec <- sprintf("%s/images/%s/input.tif", output_final, basename(s2_id))
    cloudmask_spec <- sprintf("%s/images/%s/target.tif", output_final, basename(s2_id))

    dir.create(output_final, showWarnings = FALSE)
    dir.create(output_final_d, showWarnings = FALSE)
    dir.create(output_final_folder, showWarnings = FALSE)

    # Create JSON
    ee_create_cloudseg(path = metadata_final)
    ee_create_metadata(
      id = basename(s2_id),
      point = c(jsonfile_r$y, jsonfile_r$x),
      path = metadata_spec
    )
    nlen <- length(names(final_stack))

    input_data <- raster::stack(
      final_stack[[1:13]]/10000, final_stack[[14]], final_stack[[15:17]], final_stack[[22:23]]
    )
    # 18-19 -> cmask_s2cloudness| cmask_s2cloudness_reclass (0,1)
    # 20-21 -> cmask_sen2cor | cmask_sen2cor_reclass (0,1,2)
    benchmarch_data <- final_stack[[18:21]]
    writeRaster(x = input_data, filename = inputdata_spec, overwrite = TRUE)
    writeRaster(x = benchmarch_data, filename = cloudmask_spec, overwrite = TRUE)
  }
}

ee_get_s1 <- function(point, s2_date, range = 60, exclude = NULL) {
  # 1. Defining range and ref kernel
  ee_new_kernel <- point$buffer(10*255)$bounds()
  s1_date_search <- list(
    init_date = (s2_date - lubridate::days(range)) %>% rdate_to_eedate(),
    last_date = (s2_date + lubridate::days(range)) %>% rdate_to_eedate()
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
      message("No Sentinel-1 image was found in a range of 2 months")
      stop(e)
    }
  )

  # 4. get number of pixels
  npixels_data <- s1_grd %>%
    ee$ImageCollection$map(
      function(img) {
        ee_reducer <- ee$Reducer$count()
        prop <- ee$Image$reduceRegion(
          image = img$select("VH"),
          reducer = ee_reducer,
          geometry = ee_new_kernel
        )
        img %>% ee$Image$set(list(npixels = prop))
      }
    ) %>%
    ee$ImageCollection$aggregate_array("npixels") %>%
    ee$Array$getInfo() %>%
    unlist() %>%
    as.numeric()

  # 5. Only full scene 512x512
  s1_grd_id$n_pixels <- npixels_data
  s1_grd_id_filter <- s1_grd_id %>% filter(n_pixels > 250000)
  if (!is.null(exclude)) {
    s1_grd_id_filter <- s1_grd_id_filter[!(s1_grd_id_filter$id %in% exclude),]
  }

  # 6. Get the nearest image
  row_position <- which.min(abs(s1_grd_id_filter$time_start - s2_date))
  s1_grd_id_filter[row_position,][["id"]]
}


ee_merge_s2_full <- function(s2_id, s1_id, s2_date) {
  year <- format(as.Date(s2_date), "%Y-01-01") %>% as.Date()
  if (year == as.Date("2020-01-01")) {
    year <- as.Date("2019-01-01")
  }
  year_chr <- c(year - 1 , year + 1) %>% as.character()

  # 1. Create a S2 ImageCollection and filter by space and time.
  s1_grd <- ee$Image(s1_id)
  s2_2a <- ee$Image(sprintf("COPERNICUS/S2_SR/%s", basename(s2_id)))
  s2_1c <- ee$Image(sprintf("COPERNICUS/S2/%s", basename(s2_id)))
  extra_dem <- ee$Image("MERIT/Hydro/v1_0_1")$select("elv")
  extra_LC <- ee$ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global") %>%
    ee$ImageCollection$filterDate(year_chr[1] , year_chr[2]) %>%
    ee$ImageCollection$first() %>%
    ee$Image$select("discrete_classification") %>%
    ee$Image$rename("land_cover")
  # 2. Create a S2_CLOUD_PROBABILITY ImageCollection filtering by space and time.
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

  # 3. Estimate CDI
  s2_cdi <- ee$Algorithms$Sentinel2$CDI(s2_1c)

  # 4. Estimate CDI
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

ee_create_cloudseg <- function(path) {
  cseg_list <- list(
    name = "cloud-segmentation",
    authentication_required = TRUE,
    images = list(
      path = list(
        Sentinel2 = "images/{id}/input.tif",
        CloudMask = "images/{id}/target.tif"
      ),
      shape = c(511,511),
      thumbnails = "images/{id}/thumbnail.png",
      metadata = "images/{id}/metadata.json"
    ),
    segmentation = list(
      path = "images/{id}/{id}.png",
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
        data = "$Sentinel2.B11**0.8*5",
        cmap = "jet"
      ),
      cloud_index = list(
        description = "Cloud Displacement Index, clouds are red.",
        type = "image",
        data = "$Sentinel2.B14*-1"
      ),
      "Cirrus-Edges" = list(
        "description" = "Edges in the cirrus band",
        "type" = "image",
        "data" = "edges($Sentinel2.B11**0.8*5)*1.5",
        "cmap" = "gray"
      ),
      RGB = list(
        "description" = "Normal RGB image.",
        "type" = "image",
        "data" = c("$Sentinel2.B5", "$Sentinel2.B3", "$Sentinel2.B2")
      ),
      NRGB = list(
        description = "Near-Infrared RGB image.",
        type = "image",
        data = c("$Sentinel2.B5*1.5", "$Sentinel2.B3*1.5", "$Sentinel2.B2*1.5")
      ),
      Edges = list(
        description = "Edges in the panchromatic bands",
        type = "image",
        data = "edges($Sentinel2.B2+$Sentinel2.B3+$Sentinel2.B4)",
        cmap = "gray"
      ),
      Snow = list(
        description = "Small ice crystals in high-level clouds appear reddish-orange or peach, and thick ice snow looks vivid red (or red-orange). Bare soil appears bright cyan and vegetation seem greenish in the image. Water on the ground is very dark as it absorbs the SWIR and the red, but small (liquid) water drops in the clouds scatter the light equally in both visible and the SWIR, and therefore it appears white. Water Sediments are displayed as dark red.",
        type = "image",
        data = c("$Sentinel2.B1", "$Sentinel2.B12", "$Sentinel2.B13")
      ),
      "Sentinel-1" = list(
        description = "RGB of VH, VV and VH-VV.",
        type = "image",
        data = c("$Sentinel2.B16", "$Sentinel2.B15", "$Sentinel2.B16-$Sentinel2.B15")
      ),
      Superpixels = list(
        description = "Superpixels in the panchromatic bands",
        type = "image",
        data = "superpixels($Sentinel2.B2+$Sentinel2.B3+$Sentinel2.B4, sigma=4, min_size=100)",
        cmap = "jet"
      ),
      Bing = list(
        description = "Aerial Imagery",
        type = "bingmap"
      ),
      elevation = list(
        description = "Elevation values",
        type = "image",
        data = "$Sentinel2.B18"
      ),
      sen2cloudness = list(
        description = "Sen2Cloudness Probability",
        type = "image",
        data = "$CloudMask.B1"
      ),
      sen2cloudness_reclass = list(
        description = "Sen2Cloudness Probability Reclass (BLUE-->CLEAR; RED -> CLOUD)",
        type = "image",
        data = "$CloudMask.B2"
      ),
      sen2cor = list(
        "description" = "sen2cor classes (BLUE--> CLEAR; GREEN -> CLOUD; RED -> CLOUD SHADOW)",
        "type" = "image",
        "data" = "$CloudMask.B4"
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
