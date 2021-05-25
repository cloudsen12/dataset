#' Global mean and sd from cloudSEN12
#' @author Cesar Aybar

# libraries
library(reticulate)
library(raster)
library(dplyr)
np <- import("numpy")


# Dataset data_dir
DATASET_DIR <- "/media/csaybar/58059B472A3AA231/minicloudSEN12/"
input_files <- list.files(DATASET_DIR, "input\\.npy", recursive = TRUE, full.names = TRUE)

# Estimate the mean dataset
sum_container <- rep(0, 19)
for(index in seq_along(input_files)) {
  fnp <- np$load(input_files[index])
  img_sum <- sapply(1:19, function(x) sum(fnp[x, 1:511, 1:511], na.rm = TRUE))
  sum_container <- sum_container + img_sum
}
total_mean <- sum_container/ (length(input_files)*511*511)


# Estimate the sd dataset
sum_square_container <- rep(0, 19)
for(index in seq_along(input_files)) {
  fnp <- np$load(input_files[index])
  img_sum_sq <- sapply(1:19, function(x) sum((fnp[x, 1:511, 1:511] - total_mean[x])**2, na.rm = TRUE))
  sum_square_container <- sum_square_container + img_sum_sq
}
total_sd <- sqrt(sum_square_container / (length(input_files)*511*511 - 1))


# Adding Popular Index ----------------------------------------------------

# Estimate the mean dataset
sum_container <- rep(0, 3)
for(index in seq_along(input_files)) {
  fnp <- np$load(input_files[index])
  mrange <- 1:511
  hot <- (fnp[1, mrange, mrange] - fnp[3, mrange, mrange]*0.5 - 0.08)
  hot_sum <- sum(hot, na.rm = TRUE)
  ndsi <- (fnp[2, mrange, mrange] - fnp[11, mrange, mrange])/(fnp[2, mrange, mrange] + fnp[11, mrange, mrange])
  ndsi_sum <-  sum(ndsi, na.rm = TRUE)
  ndvi <- (fnp[7, mrange, mrange] - fnp[3, mrange, mrange])/(fnp[7, mrange, mrange] + fnp[3, mrange, mrange])
  ndvi_sum <-  sum(ndvi, na.rm = TRUE)
  sum_container <- sum_container + c(hot_sum, ndsi_sum, ndvi_sum)
}
total_mean_index <- sum_container/ (length(input_files)*511*511)


# Estimate the sd dataset
sum_square_container <- rep(0, 3)
for(index in seq_along(input_files)) {
  fnp <- np$load(input_files[index])
  mrange <- 1:511
  hot <- (fnp[1, mrange, mrange] - fnp[3, mrange, mrange]*0.5 - 0.08)
  hot_sd <- sum((hot - total_mean_index[1])**2, na.rm=TRUE)

  ndsi <- (fnp[2, mrange, mrange] - fnp[11, mrange, mrange])/(fnp[2, mrange, mrange] + fnp[11, mrange, mrange])
  ndsi_sd <- sum((ndsi - total_mean_index[2])**2, na.rm=TRUE)

  ndvi <- (fnp[7, mrange, mrange] - fnp[3, mrange, mrange])/(fnp[7, mrange, mrange] + fnp[3, mrange, mrange])
  ndvi_sd <- sum((ndvi - total_mean_index[3])**2, na.rm=TRUE)

  sum_square_container <- sum_square_container + img_sum_sq
}
total_sd_index <- sqrt(sum_square_container / (length(input_files)*511*511 - 1))


# Create a table with the mean and sd values
names_in_order <- c(
  "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11",
  "B12", "VH", "VV", "angle", "CDI", "cloudshadow_direction", "elevation", "landuse",
  "hot","ndvi", "ndsi"
)

dt_sdmean <- tibble(
  name = names_in_order,
  index = 0:22,
  description = NA,
  mean = c(total_mean, NA, ),
  sd = c(total_sd, NA, )
)

# write the table
write.csv(
  x = dt_sdmean,
  file = "/home/csaybar/Documents/Github/cloudsen12/dataset/data/input_sum_table.csv",
  row.names = FALSE
)
