library(tidyverse)
library(googledrive)

## upload folder and its contents
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
