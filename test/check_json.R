library(reprex)
library(tidyverse)
library(jsonlite)

# test01 -- json format
json_control_1 <- function(metadata_files) {
  sapply(
    X = metadata_files,
    FUN = function(x) tryCatch(
      expr = {jsonlite::read_json(x); FALSE},
      error = function(e) TRUE
    )
  ) %>% as.logical()
}

#test02 -- empty JSON
json_control_2 <- function(metadata_files) {
  lapply(
    X = metadata_files,
    FUN = function(x) tryCatch(
      expr = {jsonlite::read_json(x) %>% names()},
      error = function(e) TRUE
    )
  )
}

#test03 -- Read comments
json_control_3 <- function(metadata_files) {
  lapply(
    X = metadata_files,
    FUN = function(x) tryCatch(
      expr = {
        x_com <- jsonlite::read_json(x)[["comments"]]
        if (x_com == "PUT_HERE_YOUR_COMMENT") {
          NULL
        } else {
          sprintf("%s: %s",basename(x),  x_com)
        }
        },
      error = function(e) TRUE
    )
  )
}

# duplicated ID?
json_control_4 <- function(metadata_files) {
  lapply(
    X = metadata_files,
    FUN = function(x) tryCatch(
      expr = { any(duplicated(names(jsonlite::read_json(x))))},
      error = function(e) TRUE
    )
  ) %>% unlist()
}


# The S2 tile exist?
json_control_5 <- function(points) {
  point_list <- list()
  counter <- 1
  for (point in points) {
    metadata_json <- sprintf("metadata_%04d.json", point)
    httr::set_config(httr::config( ssl_verifypeer = 0L))
    httr::set_config(httr::config(http_version = 0))
    jsonfile <- try(search_metajson(pattern = metadata_json))

    # 1. Read JSON file
    jsonfile_r <- try(jsonlite::read_json(jsonfile))
    if (class(jsonfile_r) == "try-error") {
      next
    }

    # 2. Identify all the S2 images
    s2_idsposition <- which(sapply(strsplit(names(jsonfile_r), "_"), length) == 3)
    s2_ids <- sprintf("COPERNICUS/S2/%s", names(jsonfile_r)[s2_idsposition])

    point_list[[counter]] <- data_frame(
      point = sprintf("point_%04d", point),
      image = s2_ids,
      exist = sapply(s2_ids, function(x) class(try(ee$Image(x)$getInfo())) == "try-error")
    )
    counter <- counter + 1
  }
  bind_rows(point_list)
}

#TEST ID
metadata_folder <- "metadata/"
metadata_files <- list.files(metadata_folder, full.names = TRUE)
metadata_files[json_control_1(metadata_files)]

#TEST NAME
metadata_folder <- "metadata/"
metadata_files <- list.files(metadata_folder, full.names = TRUE)
metadata_files %>%
  json_control_2 %>%
  sapply(function(x) any(grepl("PUT_HERE_ID", x))) %>%
  which() -> id_error
metadata_files[id_error]

#TEST READ_COMMENTS
metadata_folder <- "metadata/"
metadata_files <- list.files(metadata_folder, full.names = TRUE)
metadata_files %>%
  json_control_3 %>%
  unlist()

#TEST DUPLICATED ID
metadata_files[json_control_4(metadata_files)]


# Test if s2 exist
ee_Initialize("lesly")
db_points <- json_control_5(1:1500)

