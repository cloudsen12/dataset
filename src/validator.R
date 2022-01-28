library(magick)
library(raster)
library(rgee)
library(sf)

ee_Initialize()

png4bands <- function(name, alpha = 1) {
  png_mtx <- png::readPNG(name)
  if (dim(png_mtx)[3] == 3) {
    png_mtx <- abind::abind(png_mtx, matrix(alpha , 511, 511), along =3)
  }
  if (dim(png_mtx)[3] == 2) {
    png_mtx <- abind::abind(png_mtx, matrix(alpha , 511, 511), along =3)
    png_mtx <- abind::abind(png_mtx, matrix(alpha , 511, 511), along =3)
  }
  png::writePNG(png_mtx, name)
}

generate_validator_png <- function(name, ouput) {
  # Change white by transparent
  img_init <- image_read(name)
  png_tmp00 <- tempfile(fileext = ".png")
  system(sprintf("convert %s -transparent white  %s", name, png_tmp00))

  # 1. GET basename png image.
  fullname_str <- basename(name)

  # 2. Split basename you got point and image name.
  point_img_name <- strsplit(fullname_str, "__")[[1]]

  # 3. Get fullname in cloudsen12 folder
  high_img_folder <- "/media/csaybar/Elements SE/cloudSEN12/high"
  files <- sprintf("%s/%s/%s/input", high_img_folder, point_img_name[1], point_img_name[2])
  r_RGB <- stack(sprintf("%s/%s.tif", files, c("B4", "B3", "B2")))
  r_cirrus <- raster(sprintf("%s/%s.tif", files, "B10"))

  #thumbnail
  thumbpng <- stack(sprintf("%s/thumbnails/thumbnail.tif", dirname(files)))

  # 4. Save a sentinel-2 TOA png file
  png_tmp01 <- tempfile(fileext = ".png")
  png(png_tmp01, height = 511, width = 511)
  plotRGB(r_RGB, stretch = "lin")
  dev.off()

  # 5. Get the study area
  region <- as(extent(r_RGB), 'SpatialPolygons') %>% st_as_sfc()
  st_crs(region) <- as.character(crs(r_RGB))

  # 5. Save Sentinel-2 RGB BOA
  s2sr_img <- ee_as_thumbnail(
    image = ee$Image(paste0("COPERNICUS/S2_SR/", point_img_name[2]))$select(c("B4", "B3", "B2")),
    region = region %>% sf_as_ee(),
    dimensions = c(511, 511),
    vizparams = c(min = 0, max = 14000),
    raster = TRUE
  )

  png_tmp02 <- tempfile(fileext = ".png")
  png(png_tmp02, height = 511, width = 511)
  plotRGB(s2sr_img,  margins = FALSE, stretch = "hist")
  dev.off()
  #system(sprintf("convert %s -transparent white  %s", png_tmp02, png_tmp02))


  # 6. Save a sentinel-2 cirrus png file
  png_tmp03 <- tempfile(fileext = ".png")
  png(png_tmp03, height = 511, width = 511)
  graphics::par(plt=c(0,1,0,1))
  image(r_cirrus)
  dev.off()


  crs(thumbpng) <- NA
  extent(thumbpng) <- c(0,1021,0,1021)
  poly <- st_multilinestring(
    list(matrix(
      c(511/2, 511/2, 511/2, 511 + 511/2, 511 + 511/2,
        511 + 511/2, 511 + 511/2, 511/2, 511/2, 511/2),ncol = 2,byrow = T
    ))
  ) %>% as("Spatial")

  png_tmp04 <- tempfile(fileext = ".png")
  png(png_tmp04, height = 511, width = 511)
  plotRGB(thumbpng,  margins = FALSE, stretch = "lin")
  plot(poly, lwd = 2,col = "red", add = TRUE)
  dev.off()

  png_tmp05 <- tempfile(fileext = ".png")
  png(png_tmp05, height = 511, width = 511)
  plotRGB(thumbpng,  margins = FALSE, stretch = "hist")
  plot(poly, lwd = 2,col = "red", add = TRUE)
  dev.off()

  # TARGET+RGB+RGB2+cirrus
  lapply(list(png_tmp00, png_tmp01, png_tmp02, png_tmp03, png_tmp04, png_tmp05), png4bands)
  a_transparent <- image_read(png_tmp00)
  amtx <- a_transparent[[1]]
  amtx[4,,] <-  as.raw(as.integer(amtx[4,,]) * 0.3)
  amtx <- image_read(amtx)

  b <- image_read(png_tmp01)
  c <- image_read(png_tmp02)
  d <- image_read(png_tmp03)
  e <- image_mosaic(c(c, amtx))
  f <- image_read(png_tmp04)
  g <- image_read(png_tmp05)

  both <- image_scale(c(b, c, d, e, f, g), "511")
  h <- image_append(both[c(4, 3, 5)])
  i <- image_append(both[c(1, 2, 6)])
  k <- image_append(c(h, i), stack = TRUE)

  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  output_name1 <- paste0(output, "/", point_img_name[1], "/", point_img_name[2], "__01.png")
  dir.create(dirname(output_name1))
  image_write(k, output_name1)

  output_name2 <- paste0(output, "/", point_img_name[1], "/", point_img_name[2], "__02.png")
  dir.create(dirname(output_name2))
  image_write(a_transparent, output_name2)
}

# RGB + TARGET
name <- "/home/csaybar/Desktop/cloudsen12_val/progress/jhomira/"
output <- "/home/csaybar/Desktop/cloudsen12_val/validator/jhomira/"


fullname_files <- list.files(name, "\\.png$",full.names = TRUE)
for (index in 580:length(fullname_files)) {
  outputf <- paste0(output, strsplit(basename(fullname_files[index]), "__")[[1]][1])
  generate_validator_png(name = fullname_files[index], ouput = outputf)
}
