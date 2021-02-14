library(ggplot2)
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
    filter(date > as.Date("2021-02-12")) %>%
    group_by(author) %>%
    summarise(labels = length(author))
  ggplot(img_sel_df, aes(x = author, y = labels)) +
    geom_bar(stat="identity") +
    theme_classic()
}

tile_monitoring()
