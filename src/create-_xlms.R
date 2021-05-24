library(dplyr)
USER_FOLDER <- "/media/csaybar/58059B472A3AA231/angie_points/"


full_files <- list.files(USER_FOLDER,full.names = TRUE)
#full_files <- c(full_files[25:145], full_files[1:24])
full_df <- lapply(seq_along(full_files), function(x) {
  tibble(
    point = basename(full_files[x]),
    type = "high_quality",
    sen2_id = list.files(full_files[x])[1:5]
  )
})
write_csv(bind_rows(full_df), "/home/csaybar/Desktop/user_points.csv")
