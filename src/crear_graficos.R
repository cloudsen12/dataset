library(googlesheets4)
library(googledrive)
library(reticulate)
library(magick)
library(raster)
library(rgee)
library(sf)

np <- import("numpy")

tempfile <- function(case = 1, fileext = ".png") {
  sprintf("/home/csaybar/ddddx_%02d%s", case, fileext)
}

png4bands <- function(name, alpha = 1) {
  png_mtx <- png::readPNG(name)
  if (dim(png_mtx)[3] == 3) {
    png_mtx <- abind::abind(png_mtx, matrix(alpha , 509, 509), along =3)
  }
  if (dim(png_mtx)[3] == 2) {
    png_mtx <- abind::abind(png_mtx, matrix(alpha , 509, 509), along =3)
    png_mtx <- abind::abind(png_mtx, matrix(alpha , 509, 509), along =3)
  }
  png::writePNG(png_mtx, name)
}


png4bands2 <- function(name, alpha = 1) {
  png_mtx <- png::readPNG(name)
  png_mtx <- png_mtx[2:510, 2:510, 1:3]
  if (dim(png_mtx)[3] == 3) {
    png_mtx <- abind::abind(png_mtx, matrix(alpha , 509, 509), along =3)
  }
  if (dim(png_mtx)[3] == 2) {
    png_mtx <- abind::abind(png_mtx, matrix(alpha , 509, 509), along =3)
    png_mtx <- abind::abind(png_mtx, matrix(alpha , 509, 509), along =3)
  }
  png::writePNG(png_mtx, name)
}

generate_validator_png <- function(name, point, output) {
  # Change white by transparent
  img_init <- image_read(name)
  png_tmp00 <- tempfile(1, fileext = ".png")
  system(sprintf("convert %s -transparent white  %s", name, png_tmp00))
  
  # 1. GET basename png image.
  fullname_str <- basename(name)
  
  # 2. Split basename you got point and image name.
  point_img_name <- strsplit(fullname_str, "__")[[1]]
  
  # ref_raster 
  rref <- raster(list.files(point, "sen2cor\\.tif$", full.names = TRUE, recursive = TRUE)[1])
  
  # 3. Get fullname in cloudsen12 folder
  high_img_folder <- dirname(point)
  npfile <- sprintf("%s/%s/%s/input.npy", high_img_folder, point_img_name[1], point_img_name[2])
  file <- np$load(npfile)
  
  # RGB
  rref[] <- file[4,,] %>% t
  b4 <- rref/10000
  rref[] <- file[3,,] %>% t
  b3 <- rref/10000
  rref[] <- file[2,,] %>% t
  b2 <- rref/10000
  r_RGB <- stack(b4, b3, b2)
  
  # SNOW
  rref[] <- file[1,,] %>% t
  b1 <- rref/10000
  rref[] <- file[12,,] %>% t
  b11 <- rref/10000
  rref[] <- file[13,,] %>% t
  b12 <- rref/10000
  r_snow <- stack(b1, b11, b12)
  
  newpng_tmp01 <- tempfile(2, fileext = ".png")
  png(newpng_tmp01, height = 509, width = 509)
  plotRGB(r_snow,  margins = FALSE, stretch = "lin")
  dev.off()
  
  # HOT
  newpng_tmp02 <- tempfile(3, fileext = ".png")
  png(newpng_tmp02, height = 509, width = 509)
  plotRGB(r_RGB, stretch = "hist")
  dev.off()
  
  # CIRRUS
  rref[] <- file[11,,] %>% t
  r_cirrus <- rref/10000
  
  # Thumbnail
  thumbpng <- stack(
    sprintf("%s/%s/%s/thumbnail.tif", high_img_folder, point_img_name[1], point_img_name[2])
  )
  
  # 4. Save a sentinel-2 TOA png file
  png_tmp01 <- tempfile(4, fileext = ".png")
  png(png_tmp01, height = 509, width = 509)
  plotRGB(r_RGB, stretch = "lin")
  dev.off()
  
  # 5. Get the study area
  region <- as(extent(r_RGB), 'SpatialPolygons') %>% st_as_sfc()
  st_crs(region) <- as.character(crs(r_RGB))
  
  # 6. Save a sentinel-2 cirrus png file
  png_tmp03 <- tempfile(5, fileext = ".png")
  png(png_tmp03, height = 509, width = 509)
  graphics::par(plt=c(0,1,0,1))
  image(r_cirrus)
  dev.off()
  
  
  crs(thumbpng) <- NA
  extent(thumbpng) <- c(0, 951, 0, 951)
  corner <- 951/2 - 509/2
  poly <- st_multilinestring(
    list(matrix(
      c(corner, corner, corner, 509 + corner, 509 + corner,
        509 + corner, 509 + corner, corner, corner, corner),ncol = 2,byrow = T
    ))
  ) %>% as("Spatial")
  
  png_tmp04 <- tempfile(6, fileext = ".png")
  png(png_tmp04, height = 509, width = 509)
  plotRGB(thumbpng,  margins = FALSE, stretch = "lin")
  plot(poly, lwd = 2,col = "red", add = TRUE)
  dev.off()
  
  png_tmp05 <- tempfile(7, fileext = ".png")
  png(png_tmp05, height = 509, width = 509)
  plotRGB(thumbpng,  margins = FALSE, stretch = "hist")
  plot(poly, lwd = 2,col = "red", add = TRUE)
  dev.off()
  
  # TARGET+RGB+RGB2+cirrus
  lapply(list(png_tmp00, png_tmp01, png_tmp03, png_tmp04, png_tmp05, newpng_tmp01, newpng_tmp02), png4bands)
  a_transparent <- image_read(png_tmp00)
  amtx <- a_transparent[[1]]
  amtx[4,,] <-  as.raw(as.integer(amtx[4,,]) * 0.3)
  amtx <- image_read(amtx)
  
  b <- image_mosaic(c(image_read(png_tmp01), amtx))
  d <- image_read(png_tmp03)
  
  e <- image_read(png_tmp04)
  f <- image_read(png_tmp05)
  g <- image_read(newpng_tmp01)
  h <- image_mosaic(c(image_read(newpng_tmp02), amtx))
  
  both <- image_scale(c(b, d, e, f, g, h), "509")
  h <- image_append(both[c(4, 6, 5)])
  i <- image_append(both[c(3, 1, 2)])
  k <- image_append(c(h, i), stack = TRUE)
  
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  output_name1 <- paste0(output, "/", point_img_name[1], "/", point_img_name[2], "__01.png")
  dir.create(dirname(output_name1), showWarnings = FALSE, recursive = TRUE)
  image_write(k, output_name1)
  
  output_name2 <- paste0(output, "/", point_img_name[1], "/", point_img_name[2], "__02.png")
  dir.create(dirname(output_name2),showWarnings = FALSE, recursive = TRUE)
  image_write(a_transparent, output_name2)
}


all_png <- list.files(
  path = "/home/csaybar/Downloads/hq_labels_manual/jhomira/",
  pattern = "\\.png$",
  recursive = TRUE,
  full.names = TRUE
)
fixing <- lapply(all_png, png4bands2)


# (length(all_png)+600)/10000*100
point <- list.files("/media/csaybar/Elements/cloudSEN12/high/",full.names = TRUE)
output <- "/home/csaybar/Desktop/to_validate/"


# 664
points <- unique(unlist(lapply(strsplit(basename(all_png), "__"), function(x) x[1])))
for (index in 1:length(points)) {
  pname <- point[basename(point) %in% points[index]]
  name <- all_png[grepl(basename(pname), all_png)]
  if (length(name) == 0) {
    next
  }
  for (index2 in 1:5) {
    namex <- name[index2]
    generate_validator_png(name = namex, point = pname, output = output) 
  }
  print(index)
}
