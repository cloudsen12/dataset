library(ggplot2)
library(lubridate)
library(tidyverse)
library(googledrive)
library(googlesheets4)

googledrive::drive_auth("s1078735@stud.sbg.ac.at")
googlesheets4::gs4_auth("s1078735@stud.sbg.ac.at")

# monitoring progress :3
tile_monitoring <- function(week) {
  drive_jsonfile <- drive_ls(
    path = as_id("133GLIclZiNRqBfLA7aEInXyOI-sj1bYB"),
    recursive = TRUE
  )
  zips <- drive_jsonfile[grepl("\\.zip$", drive_jsonfile$name),]
  zips$name <- tolower(zips$name)
  tmpf <- tempfile(fileext = ".zip")

  ## Lissette
  tmpfolder <- sprintf("%s/images", dirname(tmpf))
  system(sprintf("rm -R %s", tmpfolder))
  dir.create(tmpfolder, showWarnings = FALSE)

  lissette <- zips[grepl("lissette", zips$name),]
  googledrive::drive_download(file = lissette, path = tmpf)
  unzip(zipfile <- tmpf, exdir = tmpfolder)
  lissette_numb <- length(list.files(tmpfolder, recursive = TRUE))/5
  lissette_fn <- lissette_numb/115*100
  message("lissette_ok")

  ## Angie
  tmpfolder <- sprintf("%s/images", dirname(tmpf))
  system(sprintf("rm -R %s", tmpfolder))
  dir.create(tmpfolder, showWarnings = FALSE)

  angie <- zips[grepl("angie", zips$name),]
  googledrive::drive_download(file = angie, path = tmpf, overwrite = TRUE)
  unzip(zipfile <- tmpf, exdir = tmpfolder)
  angie_numb <- length(list.files(tmpfolder, recursive = TRUE))/5
  angie_fn <-  angie_numb/114*100
  message("angie_ok")

  ## nicole
  tmpfolder <- sprintf("%s/images", dirname(tmpf))
  system(sprintf("rm -R %s", tmpfolder))
  dir.create(tmpfolder, showWarnings = FALSE)

  nicole <- zips[grepl("nicole", zips$name),]
  googledrive::drive_download(file = nicole, path = tmpf, overwrite = TRUE)
  unzip(zipfile <- tmpf, exdir = tmpfolder)
  nicole_numb <- length(list.files(tmpfolder, recursive = TRUE))/5
  nicole_fn <- nicole_numb/127*100
  message("nicole_ok")

  ## lesly
  tmpfolder <- sprintf("%s/images", dirname(tmpf))
  system(sprintf("rm -R %s", tmpfolder))
  dir.create(tmpfolder, showWarnings = FALSE)

  lesly <- zips[grepl("lesly|les", zips$name),]
  googledrive::drive_download(file = lesly, path = tmpf, overwrite = TRUE)
  unzip(zipfile <- tmpf, exdir = tmpfolder)
  lesly_numb <- length(list.files(tmpfolder, recursive = TRUE))/5
  lesly_fn <- lesly_numb/100*100
  message("lesly_ok")

  ## roy
  tmpfolder <- sprintf("%s/images", dirname(tmpf))
  system(sprintf("rm -R %s", tmpfolder))
  dir.create(tmpfolder, showWarnings = FALSE)

  roy <- zips[grepl("roy", zips$name),]
  googledrive::drive_download(file = roy, path = tmpf, overwrite = TRUE)
  unzip(zipfile <- tmpf, exdir = tmpfolder)
  roy_numb <- length(list.files(tmpfolder, recursive = TRUE))/5
  roy_fn <- roy_numb/90*100
  message("roy_ok")

  ## Valeria
  tmpfolder <- sprintf("%s/images", dirname(tmpf))
  system(sprintf("rm -R %s", tmpfolder))
  dir.create(tmpfolder, showWarnings = FALSE)

  val <- zips[grepl("val", zips$name),]
  googledrive::drive_download(file = val, path = tmpf, overwrite = TRUE)
  unzip(zipfile <- tmpf, exdir = tmpfolder)
  val_numb <- length(list.files(tmpfolder, recursive = TRUE))/5
  valeria_fn <- val_numb/120*100
  message("valeria_ok")

  df_progress <- data.frame(
    name = c("lissette", "angie", "nicole", "lesly", "roy", "valeria"),
    percentage = c(lissette_fn, angie_fn, nicole_fn, lesly_fn, roy_fn, valeria_fn),
    week = week
  )
}

library(readr)
dfweek <- tile_monitoring(week = 8)

# write_csv(dfweek, file = "/home/csaybar/Desktop/control.csv")
control_db <- read_csv("/home/csaybar/Desktop/control.csv")

#control_db <- control_db[control_db$name != "Prudencio",]
# control_db <- control_db[control_db$name == "angie",]
control_db <- rbind(control_db, dfweek)
# write_csv(control_db, file = "/home/csaybar/Desktop/control.csv")

ggplot(control_db, aes(x = factor(week), y = percentage, group = name, color = name)) +
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  ylim(c(0, 110)) +
  geom_hline(yintercept = 100, color = "red") +
  xlab("WEEK")
