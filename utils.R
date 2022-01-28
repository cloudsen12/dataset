# Land cover classes (LC10)
LAND_COVER_CLASSES <- list(
  "10" = "Trees",
  "20" = "Shrubland",
  "30" = "Grassland",
  "40" = "Cropland",
  "50" = "Built-up",
  "60" = "Barren / sparse vegetation",
  "70" = "Snow and ice",
  "80" = "Open water",
  "90" = "Herbaceous wetland",
  "95" = "Mangroves",
  "100" = "Moss and lichen"
)

ip_creator <- function(dataset) {
  # Directory path
  dirpath <- sprintf("%s/%s", CLOUDSEN12_PATH, dataset$label_type)

  # Load S2 images
  s2l1c <- ee$Image(sprintf("COPERNICUS/S2/%s", dataset$sen2))
  s2l2a <- ee$Image(sprintf("COPERNICUS/S2_SR/%s", dataset$sen2))
  s2_date <- ee_get_date_img(s2l1c)$time_start
  
  # Load centroid 
  st_point_geo <- st_as_sfc(dataset$proj_centroid, crs = 4326)
  crs_kernel <- s2l1c$select(0)$projection()$getInfo()$crs
  st_point_utm <- st_transform(st_point_geo, crs_kernel)
  ee_point <- ee$Geometry$Point(st_point_utm[[1]], proj = crs_kernel)
  
  # Get the closest S1 image.
  s1_id <- ee_get_s1(point = ee_point, s2_date = s2_date, range = 2.5)
  s1 <- ee$Image(s1_id)
  
  # extra
  s2_cdi <- ee$Algorithms$Sentinel2$CDI(s2l1c) %>%
    ee$Image$unmask(-99, sameFootprint = F)
  s2_shadowdir <- shadow_direction(s2l1c) %>% 
    ee$Image$rename("shadowdir")
  elevation <- cloudsen12_dem() %>% 
    ee$Image$rename("elevation")
  ocurrence <- ee$Image("JRC/GSW1_3/GlobalSurfaceWater") %>% 
    ee$Image$select("occurrence") %>% 
    ee$Image$rename("jrc_water_ocurrence") %>% 
    ee$Image$unmask(-99, sameFootprint = FALSE)
  lc100 <- ee$ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global") %>% 
    ee$ImageCollection$filterDate("2019-01-01", "2019-12-31") %>% 
    ee$ImageCollection$first() %>% 
    ee$Image$select("discrete_classification") %>% 
    ee$Image$rename("landcover-100") %>% 
    ee$Image$unmask(-99, sameFootprint = FALSE)
  lc10 <- ee$ImageCollection("ESA/WorldCover/v100") %>% 
    ee$ImageCollection$filterDate("2020-01-01", "2020-12-31") %>% 
    ee$ImageCollection$first() %>% 
    ee$Image$rename("landcover-10") %>% 
    ee$Image$unmask(-99, sameFootprint = FALSE)
  extra <- ee$Image$cat(list(s2_cdi, s2_shadowdir, elevation, ocurrence, lc100, lc10))
  
  # Download extra
  message("Downloading extra/ features")
  stars_extra <- fast_from_ee_to_local(
    image = extra,
    crs_kernel = crs_kernel,
    ee_point = ee_point
  )
  
  # Download s2l1c
  message("Downloading S2L1C features")
  s2l1c_extra <- fast_from_ee_to_local(
    image = s2l1c,
    crs_kernel = crs_kernel,
    ee_point = ee_point
  )
  
  # Download s2l2a
  message("Downloading S2L2A features")
  s2l2a_extra <- fast_from_ee_to_local(
    image = s2l1c,
    crs_kernel = crs_kernel,
    ee_point = ee_point
  )
  
  # Download s1
  message("Downloading S1 features")
  s1_extra <- fast_from_ee_to_local(
    image = s1,
    crs_kernel = crs_kernel,
    ee_point = ee_point
  )
  
  # Export final results.
  list(
    extra = stars_extra,
    s2l1c = s2l1c_extra,
    s2l2a = s2l2a_extra,
    s1 = s1_extra
  )
}

metadata_creator <- function(dataset, raster_ref) {
  # Label type
  label_type <- dataset$label_type
  
  # Directory path
  dirpath <- sprintf("%s/%s", CLOUDSEN12_PATH, dataset$label_type)
  
  # point id
  roi_id <- dataset$ROI
  
  # s2 id GEE
  sentinel2_product_id_gee <- dataset$sen2
  
  # s2 id
  ee_s2 <- ee$Image(sprintf("COPERNICUS/S2/%s", dataset$sen2))
  sentinel2_product_id <- ee_s2$get("PRODUCT_ID")$getInfo()
  
  # s2 date
  sentinel2_date <- ee_get_date_img(ee_s2)[["time_start"]]
  
  # PROJ extention ----------------------------------------------------------
  proj_geometry <- raster_ref$s2l1c %>% st_bbox() %>% st_as_sfc() %>% st_as_text()
  proj_epsg <- st_crs(raster_ref$s2l1c)$epsg
  proj_centroid <- raster_ref$s2l1c %>% st_bbox() %>% 
    st_as_sfc() %>% 
    st_centroid() %>%
    st_transform(4326) %>% 
    st_as_text()
  proj_shape <- 509
  proj_transform <- sprintf(
    "10, 0, %s, 0, -10, %s", 
    attr(raster_ref$s2l1c, "dimensions")$x$offset,
    attr(raster_ref$s2l1c, "dimensions")$y$offset
  )
  # -------------------------------------------------------------------------
  
  # s1 id
  ee_point_center <- st_as_sfc(proj_centroid) %>% st_set_crs(4326) %>% sf_as_ee()
  s1_id_gee <- ee_get_s1(ee_point_center, s2_date = sentinel2_date)
  sentinel1_product_id <- basename(s1_id_gee)
  sentinel1_date <- s1_id_gee %>% ee$Image() %>% ee_get_date_img() %>% '[['("time_start")
  
  # annotator name
  annotator_name <- dataset$user
  
  # grd_post_processing_software_name
  grd_post_processing_software_name <- s1_id_gee %>% 
    ee$Image() %>% 
    ee$Image$get("GRD_Post_Processing_software_name") %>% 
    ee$ComputedObject$getInfo()
  
  # grd_post_processing_software_version
  grd_post_processing_software_version <- s1_id_gee %>% 
    ee$Image() %>% 
    ee$Image$get("GRD_Post_Processing_software_version") %>% 
    ee$ComputedObject$getInfo()
  
  # slc_processing_facility_name
  slc_processing_facility_name <- s1_id_gee %>% 
    ee$Image() %>% 
    ee$Image$get("SLC_Processing_facility_name") %>% 
    ee$ComputedObject$getInfo()
  
  # SLC_Processing_software_version
  slc_processing_software_version <- s1_id_gee %>% 
    ee$Image() %>% 
    ee$Image$get("SLC_Processing_software_version") %>% 
    ee$ComputedObject$getInfo()
  
  # fmask_version
  fmask_version <- "4.3.0"
  
  # sencloudness_version
  sencloudness_version <- "1.5.0"
  
  # sen2cor_version
  ee_s2sr <- ee$Image(sprintf("COPERNICUS/S2_SR/%s", sentinel2_product_id_gee))
  sen2cor_dtidentf <- ee_s2sr$get("DATATAKE_IDENTIFIER")$getInfo()
  sen2cor_version <- strsplit(sen2cor_dtidentf, "_")[[1]][4]

  # VIEW extension ----------------------------------------------------------
  # view_off_nadir
  view_off_nadir <- 0
  
  # view_sun_azimuth
  view_sun_azimuth <- ee_s2 %>% 
    ee$Image$get("MEAN_SOLAR_AZIMUTH_ANGLE") %>% 
    ee$ComputedObject$getInfo()
  
  # view_sun_elevation
  view_sun_elevation <- 90 - ee_s2 %>% 
    ee$Image$get("MEAN_SOLAR_ZENITH_ANGLE") %>% 
    ee$ComputedObject$getInfo()
  # -------------------------------------------------------------------------
  
  # s1 coverage
  sar_ref <- raster_ref$s1$VV
  radar_coverage <- 1 - sum(sar_ref[[1]] == -99)/259081
  
  # land cover
  lcover <- raster_ref$extra["landcover.10"]
  lcmode <- lcover %>% '[['("landcover.10") %>% table()
  land_cover_code <- names(lcmode[which.max(lcmode)])
  land_cover_name <- LAND_COVER_CLASSES[[land_cover_code]]
  
  # reflectance_conversion_correction
  reflectance_conversion_correction <- ee_s2 %>% 
    ee$Image$get("REFLECTANCE_CONVERSION_CORRECTION") %>% 
    ee$ComputedObject$getInfo()
  
  # aot_retrieval_accuracy
  aot_retrieval_accuracy <- ee_s2sr %>% 
    ee$Image$get("AOT_RETRIEVAL_ACCURACY") %>% 
    ee$ComputedObject$getInfo()
  
  # water_vapour_retrieval_accuracy 
  water_vapour_retrieval_accuracy <- ee_s2sr %>% 
    ee$Image$get("WATER_VAPOUR_RETRIEVAL_ACCURACY") %>% 
    ee$ComputedObject$getInfo()
  
  # Wrap everything together
  tibble(
    roi_id = roi_id,
    s2_id_gee = sentinel2_product_id_gee,
    s2_id = sentinel2_product_id,
    s2_date = sentinel2_date,
    s2_sen2cor_version = sen2cor_version,
    s2_fmask_version = fmask_version,
    s2_cloudless_version = sencloudness_version,
    s2_reflectance_conversion_correction = reflectance_conversion_correction,
    s2_aot_retrieval_accuracy = aot_retrieval_accuracy,
    s2_water_vapour_retrieval_accuracy = water_vapour_retrieval_accuracy,
    s2_view_off_nadir = view_off_nadir,
    s2_view_sun_azimuth = view_sun_azimuth,
    s2_view_sun_elevation = view_sun_elevation,
    s1_id = sentinel1_product_id,
    s1_date = sentinel1_date,
    s1_grd_post_processing_software_name = grd_post_processing_software_name,
    s1_grd_post_processing_software_version = grd_post_processing_software_version,
    s1_slc_processing_facility_name = slc_processing_facility_name,
    s1_slc_processing_software_version = slc_processing_software_version,
    proj_epsg = proj_epsg,
    proj_geometry = proj_geometry,
    proj_shape = proj_shape,
    proj_centroid = proj_centroid,
    proj_transform = proj_transform,
    label_type = label_type
  )
}


#' Convert an ee.Image to an stars object
#' @param image An ee.Image
#' @param crs_kernel The coordinate reference system of the exported image's 
#' projection. Defaults to the image's default projection.
#' @param ee_point Centroid of the image patch.
fast_from_ee_to_local <- function(image, crs_kernel, ee_point) {
  # 4.7 Create a 509x509 tile (list -> tibble -> stars)
  band_names_s2 <- image$bandNames()$getInfo()
  band_names <- c(band_names_s2, "x", "y")
  s2_img_array <- image %>%
    ee$Image$addBands(ee$Image$pixelCoordinates(projection = crs_kernel)) %>%
    ee$Image$neighborhoodToArray(
      kernel = ee$Kernel$rectangle(254, 254, "pixels")
    ) %>%
    ee$Image$sampleRegions(
      collection = ee_point,
      projection = crs_kernel,
      scale = 10) %>%
    ee$FeatureCollection$getInfo()
  length(s2_img_array$features)
  
  extract_fn <- function(x) as.numeric(unlist(s2_img_array$features[[1]]$properties[x]))
  image_as_df <- do.call(cbind,lapply(band_names, extract_fn))
  colnames(image_as_df) <- band_names
  image_as_tibble <- as_tibble(image_as_df)
  as_stars <- lapply(
    X = band_names_s2, 
    FUN = function(z) st_as_stars(image_as_tibble[c("x", "y", z)])
  ) %>% do.call(c, .)
  st_crs(as_stars) <- crs_kernel
  as_stars
}


#' Get the closest Sentinel-1 image to an specific Sentinel-2 image
#' @param point The centroid of the Sentinel-2 image.
#' @param s2_date The date of the Sentinel-2 image.
#' @param range The time range in which the Sentinel-1 images are searched. Defaults to 2.5.
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
      stop(e)
    }
  )
  
  # 4. Get the nearest image
  row_position <- which.min(abs(s1_grd_id$time_start - s2_date))
  s1_grd_id[row_position,][["id"]]
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


#' Merge MERIT + CrySat2 DEM
#' @noRd
cloudsen12_dem <- function() {
  merit_dem <- ee$Image("MERIT/Hydro/v1_0_1")$select("elv")
  cryosat2_dem <- ee$Image("CPOM/CryoSat2/ANTARCTICA_DEM")$select("elevation")
  cloudsen12_dem <- ee$Image$add(merit_dem$unmask(0), cryosat2_dem$unmask(0))
  cloudsen12_dem
}
