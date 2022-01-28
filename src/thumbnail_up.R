library(doMC)
library(stars)
library(foreach)

registerDoMC(5)

thumbnails_creator <- function(point) {
  # folder1
  pointf <- sprintf("%s/high/%s", "/media/csaybar/Elements SE/cloudSEN12", point)
  ips <- list.files(pointf, full.names = TRUE) 
  to_copy <- sprintf("%s/thumbnails/thumbnail.tif", ips[dir.exists(ips)]) %>% sort()
  
  pointf2 <- sprintf("%s/high/%s", "/media/csaybar/Elements SE/cloudSEN12_f", point)
  ips2 <- list.files(pointf2, full.names = TRUE)
  to_del <-sprintf("%s/thumbnail.tif", ips2) %>% sort()
  file.copy(to_copy, to_del, overwrite = TRUE)  
  
  # crop to 1001
  allfiles <- list.files(pointf2, pattern = "\\.tif", full.names = TRUE, recursive = TRUE) 
  ref <- allfiles[grepl("sen2cor", allfiles)]
  ddd <- read_stars(ref, proxy = TRUE) # check crs
  crs <- st_crs(ddd)$epsg
  
  r <- raster(ref[1]) 
  centroid <- st_sfc(st_point(c(mean(xFromCol(r)), mean(yFromRow(r)))), crs = crs(r)@projargs)
  ref_b <- st_buffer(centroid, 500*10, endCapStyle="SQUARE")
  xmin <- as.numeric(st_bbox(ref_b))[1] - 5
  ymin <- as.numeric(st_bbox(ref_b))[2] - 5
  xmax <- as.numeric(st_bbox(ref_b))[3] + 5
  ymax <- as.numeric(st_bbox(ref_b))[4] + 5
  
  tempfiletifs <- sprintf("%s.tif", paste0("/home/csaybar/t_", sprintf("%02d",1:5)))
  foreach (index2 = 1:5) %dopar% {
    system(
      sprintf(
        "gdalwarp '%s' '%s' -overwrite -te %s %s %s %s -tr 10 10 -q -t_srs %s",
        to_del[index2], tempfiletifs[index2],
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
  read_stars(to_del, proxy = TRUE)
}

points <- list.files(sprintf("%s/high", "/media/csaybar/Elements SE/cloudSEN12_f"))
for (index in 1688:2000) {
  thumbnails_creator(point = points[index])
  print(index)
}

