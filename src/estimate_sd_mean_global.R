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


# Create a table with the mean and sd values
names_in_order <- c(
  "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11",
  "B12", "VH", "VV", "angle", "CDI", "cloudshadow_direction", "elevation", "landuse"
)

dt_sdmean <- tibble(
  name = names_in_order,
  index = 0:19,
  description = NA,
  mean = c(total_mean, NA),
  sd = c(total_sd, NA)
)

# write the table
write.csv(
  x = dt_sdmean, 
  file = "/home/csaybar/Documents/Github/cloudsen12/dataset/data/input_sum_table.csv", 
  row.names = FALSE
)
