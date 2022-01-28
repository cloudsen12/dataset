library(foreach)
library(raster)
library(stars)
library(rgee)
library(doMC)

registerDoMC(5)
ee_Initialize("gabriela")

fmask_up <- function(point) {
  fmask_folder <- "/media/csaybar/58059B472A3AA231/FMASK4"  
  allfmaskfiles <-  list.files(fmask_folder, "\\.tif$", recursive = TRUE, full.names = TRUE)
  allfmaskfiles_b <- basename(allfmaskfiles)
  fdir <- sprintf("/media/csaybar/Elements SE/cloudSEN12_f/high/%s", point)
  fimg <- list.files(fdir)
  fimg_gr_id <- sapply(1:5, 
    function(x) {
      id <- ee$Image(sprintf("COPERNICUS/S2/%s", fimg[x]))$get("GRANULE_ID")$getInfo()
      allfmaskfiles[allfmaskfiles_b %in% sprintf("%s_Fmask4.tif", id)][1]
    }
  )
  
  if (length(fimg_gr_id) != 5) {
    stop("no fmask4 img")
  }
  to_del <- sprintf("%s/%s/models/fmask.tif", fdir, fimg)
  # get ref
  # crop to 1001
  allfiles <- list.files(fdir, pattern = "\\.tif", full.names = TRUE, recursive = TRUE) 
  ref <- allfiles[grepl("sen2cor", allfiles)]
  ddd <- read_stars(ref, proxy = TRUE) # check crs
  crs <- st_crs(ddd)$epsg
  
  r <- raster(ref[1])
  rextent <- extent(r)
  xmin <- rextent[1]
  xmax <- rextent[2]
  ymin <- rextent[3]
  ymax <- rextent[4]
  
  tempfiletifs <- sprintf("%s.tif", paste0("/home/csaybar/z_", sprintf("%02d",1:5)))
  foreach (index2 = 1:5) %dopar% {
    system(
      sprintf(
        "gdalwarp '%s' '%s' -overwrite -te %s %s %s %s -tr 10 10 -q -s_srs %s",
        fimg_gr_id[index2], tempfiletifs[index2],
        xmin, ymin,
        xmax, ymax, sprintf("EPSG:%s", crs)
      )
    )
    file.copy(
      from = tempfiletifs[index2],
      to = to_del[index2],
      overwrite = TRUE
    )
  }
}

allp <- list.files("/media/csaybar/Elements SE/cloudSEN12_f/high/")
for (index in 1:2000) {
  print(index)
  fmask_up(point = allp[index]) 
}
