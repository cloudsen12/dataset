library(ggplot2)
library(lubridate)
library(tidyverse)
library(googledrive)
library(googlesheets4)

googledrive::drive_auth("s1078735@stud.sbg.ac.at")
googlesheets4::gs4_auth("s1078735@stud.sbg.ac.at")

# monitoring progress :3
tile_monitoring <- function() {
  drive_jsonfile <- drive_ls(
    path = as_id("1fBGAjZkjPEpPr0p7c-LtJmfbLq3s87RK")
  )
  see_changes_get_time <- function(x) drive_jsonfile[x,]$drive_resource[[1]]$createdTime
  see_changes_get_author <- function(x) drive_jsonfile[x,]$drive_resource[[1]]$lastModifyingUser$displayName
  author_img_selection <- sapply(1:nrow(drive_jsonfile), see_changes_get_author)
  time_img_selection <- sapply(1:nrow(drive_jsonfile), see_changes_get_time)
  remove_cesar <- !author_img_selection == "Cesar Luis Aybar"
  author_img_selection <- author_img_selection[remove_cesar] %>%
    strsplit(" ") %>%
    sapply(function(x) x[[1]][[1]])
  time_img_selection <- time_img_selection[remove_cesar] %>% as.POSIXct %>% as.character
  img_sel_df <- data_frame(author = author_img_selection, date = as.Date(time_img_selection)) %>%
    group_by(author, cut(date, "week")) %>%
    summarise(labels = length(author)) %>%
    ungroup()
  names(img_sel_df) <- c("author", "date", "nlabels")
  ggplot(img_sel_df, aes(x = date, y = nlabels, group = author, color = author)) +
      geom_line(size = 1.5) +
      geom_point(size = 2.5) +
      ylim(c(0, 300)) +
      geom_hline(yintercept = 250, color = "red")
}

# hqlabels monitoring in progress :3
hqlabels_monitoring <- function() {
  xls <- read_sheet(
    ss = as_sheets_id("1LpW9JY2BdhlQvAObD1BCzoBiliNRnU3fCRWMQJFRvoM"),
    range = "A1:F7000"
  )
  xls_db <- xls %>%
    filter(labeler  %in%  c("Jhomira", "Eduardo", "Fernando")) %>%
    filter(!is.na(sen2_id))
  txls <- table(xls_db[["labeler"]])
  df_f <- data_frame(
    labeler = names(txls),
    points = floor(txls/5)
  )

  ggplot(df_f, aes(x = labeler, y = points)) +
    geom_bar(stat="identity") +
    theme_bw()
}

tile_monitoring()
httr::set_config(httr::config( ssl_verifypeer = 0L))
httr::set_config(httr::config(http_version = 0))
hqlabels_monitoring()
