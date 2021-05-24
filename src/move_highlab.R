library(googledrive)
library(raster)
library(stars)

DIRNAME <- "/media/csaybar/Elements SE/cloudSEN12/extra/manual/"
CLOUSEN12DIR <- "/media/csaybar/Elements SE/cloudSEN12_f/high/"
CLOUSEN12DIR2 <- "/media/csaybar/Elements SE/cloudSEN12/high/"

#' From GeoTIFF to COG
create_final_target <- function(x, output) {
  clear <- x == 0
  thick_cloud <- x == 1
  thin_cloud <- x == 2
  shadow_cloud <- x == 3
  stk_stars <- stack(clear, thick_cloud, thin_cloud, shadow_cloud) %>%
    st_as_stars() %>%
    '[['("layer.1") %>%
    np$moveaxis(2L, 0L) %>%
    np$save(output, .)
  output
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

# 1. Download drive zip

## Parent folder in root folder
files <- drive_find(
  q = "'133GLIclZiNRqBfLA7aEInXyOI-sj1bYB' in parents"
)

## search for zip files and download
zip_files <- lapply(1:11, function(x) {
  zipfile <- drive_ls(path = as_id(files$id[x]))
  if (nrow(zipfile) == 0) {
    return(NA)
  }
  tmpf <- tempfile(
    tmpdir = DIRNAME,
    fileext = ".zip"
  )
  drive_download(zipfile, tmpf)
  tmpf
}) %>% unlist() %>% na.omit() %>% as.character()

# 2. Uncrompress zip files
lapply(seq_along(zip_files), function(x) unzip(zip_files[x], exdir = DIRNAME))

## Do all the points have 5 manual target png?
all_target <- list.files(DIRNAME, "\\.png$",full.names = TRUE) # all labelings
table_summ <- lapply(strsplit(basename(all_target), "__"), function(x) x[[1]]) %>%
  unlist() %>%
  table()
all(!table_summ == 5) # is it is FALSE everything going well

# 3. From extra folder to cloudSEN12

for (index in 12:400) {
  ## List file in the point folder
  llf <- list.files(CLOUSEN12DIR, full.names = TRUE)[index]
  message("Processing: ", basename(llf))
  ## Select manual labeling in point folder
  target_llf <- all_target[grepl(basename(llf), all_target)]

  ## 5 files?
  if (length(target_llf) != 5) {
    message(sprintf("%s: no 5 target images", basename(llf)))
    next
  }

  ## From extra folder to cloudsen12 folder
  for (index2 in 1:5) {
    # Get s2_id
    s2_id <- strsplit(target_llf, "__")[[index2]][2]
    ref_raster <- raster(sprintf("%s/%s/models/sen2cor.tif", llf, s2_id))
    high_labeling_r <- number_to_hex(
      stack_color = target_llf[index2],
      raster_ref = ref_raster
    )
    writeRaster(
      x = high_labeling_r,
      filename = sprintf("%s/%s/models/manual.tif", llf, s2_id),
      overwrite = TRUE
    )
    create_final_target(
      x = high_labeling_r,
      output = sprintf("%s/%s/manual.npy", llf, s2_id)
    )
  }
}


# 4. Search missing points
table_summ2 <- lapply(strsplit(basename(all_target), "__"), function(x) x[[1]]) %>%
  unlist() %>%
  unique()

CLOUSEN12DIR3 <- "/media/csaybar/Elements SE/minicloudSEN12/"
delete_points <- list.files(CLOUSEN12DIR3)[!list.files(CLOUSEN12DIR3) %in%  table_summ2]
lapply(delete_points, function(x) system(sprintf("rm -R /media/csaybar/Elements\\ SE/minicloudSEN12/%s", x)))

