library(magrittr)
library(jsonlite)

json_control <- function(metadata_files) {
  sapply(
    X = metadata_files,
    FUN = function(x) tryCatch(
      expr = {jsonlite::read_json(x); FALSE},
      error = function(e) TRUE
    )
  ) %>% as.logical()
}

metadata_folder <- "metadata/"
metadata_files <- list.files(metadata_folder, full.names = TRUE)
metadata_files[json_control(metadata_files)]

