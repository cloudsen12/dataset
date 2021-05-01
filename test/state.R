library(sf)
library(dplyr)
library(ggplot2)
library(googledrive)


# 1. functions ------------------------------------------------------------
gg_plot_ts <- function(serie, id_folder, output = "/home/csaybar/Desktop/report/") {
  ID_FOLDER_TRUE <- as.numeric(gsub("(\\D+)", "", names(id_folder[id_folder == TRUE])))
  df <- tibble(
    id = ID_FOLDER_TRUE,
    state = ID_FOLDER_TRUE %in% serie,
    number = 1
  )
  gg1 <- ggplot(df, aes(x = id, y = number, fill = state)) +
    geom_tile() +
    theme_void() +
    theme(legend.position = "none") +
    theme(
      axis.text.x = element_text(color="black", size=8, angle=90)
    ) +
    scale_x_continuous(
      limits = c(min(ID_FOLDER_TRUE), max(ID_FOLDER_TRUE)),
      breaks = seq(min(ID_FOLDER_TRUE), max(ID_FOLDER_TRUE), 100)
    )
  ggsave(sprintf("%s/ggmap_%s.png", output, max(ID_FOLDER_TRUE)), gg1)
  print(gg1)
  bad_points <- paste0(df %>% filter(state !=TRUE) %>% '[['("id"), collapse = " ")
  write.table(bad_points, sprintf("%s/bpoints_%s.txt", output, max(ID_FOLDER_TRUE)))
}

rigorous_search <- function(names) {
  save <- rep(NA, length(names))
  counter <- 1
  for (name in names) {
    folder <- googledrive::drive_ls(
      path = as_id("1NOxjbbtiyz2UqiJAJz6IsCNaLxOJlErr"),
      q = sprintf("name contains '%s'", name)
    )
    save[counter] <- nrow(folder) == 1
    counter <- counter + 1
    message(name)
  }
  names(save) <- names
  save
}

control_progress_creator <- function() {
  progress_list <- list(
    id = 1:nrow(local_cloudsen2_points),
    ts_folders = rep(NA, nrow(local_cloudsen2_points)),
    ts_json = rep(NA, nrow(local_cloudsen2_points)),
    roy_folder = rep(NA, nrow(local_cloudsen2_points)),
    final_dataset = rep(NA, nrow(local_cloudsen2_points))
  )
  write.csv(data.frame(progress_list), "data/data_process.csv", row.names = FALSE)
}


control_progress_update <- function(values, upgrade = "ts_folders") {
  if (is.null(names(values))) {
    stop("values should have a name. (e.g. point_1945)")
  }
  progress_df <- read.csv("data/data_process.csv")
  write.csv(progress_df, "data/data_process_backup.csv", row.names = FALSE)

  ID_FOLDER_TRUE <- as.numeric(gsub("(\\D+)", "", names(values)))
  progress_df[[upgrade]][ID_FOLDER_TRUE] <- values
  message(sprintf("%s column up!.", upgrade))
  write.csv(progress_df, "data/data_process.csv", row.names = FALSE)
}

# -------------------------------------------------------------------------
# control_progress_creator()

# sf dataset
local_cloudsen2_points <- read_sf("data/cloudsen2_potential_points.geojson")

# List metadata (.jsonfiles)
path <- "/home/csaybar/Desktop/dataset/metadata/"
files <- list.files(path, "\\.json$", full.names = TRUE)
metadata_id <- sort(as.numeric(gsub("metadata_|\\.json$", "", basename(files))))


for (index in 10:11) {
  # List drive folder
  is_available <- sprintf("point_%04d", (index*1000 + 1):((index + 1)* 1000))
  tile_selection_folders <- rigorous_search(is_available)

  # Plot your dataset
  gg_plot_ts(
    serie = metadata_id,
    id_folder = tile_selection_folders
  )

  # upgrade control_progress: ts_folders, ts_json, roy_folder, final_dataset
  control_progress_update(
    values = tile_selection_folders,
    upgrade = "ts_folders"
  )
}

# metadata_id <- sort(as.numeric(gsub("metadata_|\\.json$", "", basename(files))))
# values <- 1:nrow(local_cloudsen2_points) %in% metadata_id
# names(values) <- sprintf("point_%04d", 1:nrow(local_cloudsen2_points))

# control_progress_update(
#   values = values[1:1450],
#   upgrade = "ts_folders"
# )
# control_progress_update(
#   values = values,
#   upgrade = "ts_json"
# )


library(tmap)

progress_df <- read.csv("data/data_process.csv")
cloudsen12_geo <- local_cloudsen2_points[["geometry"]]
cloudsen12_geo <- st_sf(
  points = as.numeric(progress_df$ts_json),
  label = local_cloudsen2_points$label,
  size = 1,
  geometry = cloudsen12_geo
)


data("World")
data(World, metro, rivers, land)
tmap_mode("plot")

## tmap mode set to plotting
tm_shape(land) +
  tm_raster("elevation", palette = terrain.colors(10)) +
  tm_shape(World) +
  tm_borders("white", lwd = .5) +
  tm_shape(geo %>% filter(label == 'high_quality')) +
  tm_symbols(col = "points", scale = .2, border.lwd = 0,palette = c("red", "blue")) +
  tm_legend(show = FALSE)


# geo <- cloudsen12_geo %>% filter(points == 1)
# geo %>% filter(label == 'high_quality')

