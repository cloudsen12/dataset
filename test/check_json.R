library(reprex)
library(jsonlite)
library(tidyverse)
library(googledrive)

# test01 -- JSON format --- we can read the JSON file?
json_control_1 <- function(metadata_files) {
  sapply(
    X = metadata_files,
    FUN = function(x) tryCatch(
      expr = {jsonlite::read_json(x); FALSE},
      error = function(e) TRUE
    )
  ) %>% as.logical()
}

#test02 -- empty JSON --- all the JSON files have a name?
json_control_2 <- function(metadata_files) {
  lapply(
    X = metadata_files,
    FUN = function(x) tryCatch(
      expr = {jsonlite::read_json(x) %>% names()},
      error = function(e) TRUE
    )
  )
}

# test03 -- Read all the comments ---
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

# test04 -- is there duplicated ID in the same JSON?
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
json_control_5 <- function(json_name) {
    # 1. Read JSON file
    jsonfile_r <- try(jsonlite::read_json(json_name))
    if (class(jsonfile_r) == "try-error") {
      FALSE
    } else {
      # 2. Identify all the S2 images
      s2_idsposition <- which(sapply(strsplit(names(jsonfile_r), "_"), length) == 3)
      s2_ids <- sprintf("COPERNICUS/S2/%s", names(jsonfile_r)[s2_idsposition])
      any(sapply(s2_ids, function(x) class(try(ee$Image(x)$getInfo())) == "try-error"))
    }
}


full_test <- function(json_name = "metadata/metadata_0001.json") {
  # names(jsonlite::read_json(json_name))
  test_01 <- json_control_1(json_name)
  if (test_01 == TRUE) {
    test_02 <- FALSE
    # test_03 <- FALSE
    test_04 <- FALSE
    test_05 <- FALSE
  } else {
    test_02 <- any(grepl("PUT_HERE_ID", unlist(json_control_2(json_name))))
    # test_03 <- !is.null(unlist(json_control_3(json_name)))
    test_04 <- json_control_4(json_name)
    test_05 <- json_control_5(json_name)
  }
  drive_jsonfile <- drive_ls(
    path = as_id("1fBGAjZkjPEpPr0p7c-LtJmfbLq3s87RK"),
    q = sprintf("name contains '%s'", basename(json_name))
  )
  author <- drive_jsonfile$drive_resource[[1]]$lastModifyingUser$displayName
  tibble(
    json_name = basename(json_name),
    labeler = author,
    test_01 = test_01,
    test_02 = test_02,
    # test_03 = test_03,
    test_04 = test_04,
    test_05 = test_05,
  )
}



# TEST 01: Does the JSON malformed?
# TEST 02: Does the JSON have empty names?
# TEST 04: Is there duplicated Sentinel2 ID in the same JSON?
# TEST 05: Does the Sentinel2 ID malformed?

#TEST ID
setwd("/home/csaybar/Documents/Github/cloudsen12/dataset/")
metadata_folder <- "metadata/"
metadata_files <- list.files(metadata_folder, full.names = TRUE)
index <- 874
testing_json <- list()

for (index in 1:length(metadata_files)) {
  print(index)
  httr::set_config(httr::config( ssl_verifypeer = 0L))
  httr::set_config(httr::config(http_version = 0))
  testing_json[[index]] <- try(full_test(metadata_files[[index]]))
}

total_db <- bind_rows(testing_json[which(!sapply(testing_json, function(x) class(x)[[1]] == "try-error"))])
write_csv(
  x = total_db %>% filter(test_01 == TRUE | test_02 == TRUE | test_04 == TRUE | test_05 == TRUE),
  file =  "/home/csaybar/Desktop/db.csv"
)


