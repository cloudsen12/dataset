#' Rescue points!
#'
#' Sometimes the download fall :c . This script permit you merge
#' the rescue dataset with the cloudSEN12 dataset
DIRNAME <- "/home/csaybar/Desktop/cloudsen12/dataset/"


# 1. clean IPL_cloud
point_list <- list.files(
  path = DIRNAME,
  pattern = "IPL_cloudmask_reclass\\.tif$",
  full.names = TRUE,
  recursive = TRUE
)
unlink(point_list)

# 2. List points (remove points with no full data)
point_list <- list.files(DIRNAME, full.names = TRUE)
cond <- sapply(seq_along(point_list), function(x) length(list.files(point_list[x], recursive = TRUE)) == 143)
point_list[!cond]



# Parent folder in root folder
files <- drive_find(
  q = "'133GLIclZiNRqBfLA7aEInXyOI-sj1bYB' in parents"
)

# search for zip files and download
zip_files <- lapply(1:11, function(x) {
  zipfile <- drive_ls(path = as_id(files$id[x]))
  if (nrow(zipfile) == 0) {
    return(NA)
  }
  tmpf <- tempfile(
    tmpdir = "/media/csaybar/Elements SE/cloudSEN12/extra/manual/",
    fileext = ".zip"
  )
  drive_download(zipfile, tmpf)
  tmpf
}) %>% unlist() %>% na.omit() %>% as.character()
