library(jsonlite)

scribble_creator <- function(folder) {
  points <- list.files(folder, full.names = TRUE)
  message(sprintf("%s fueron detectados", length(points)))
  Sys.sleep(2)
  for (index in seq_along(points)) {
    in_files <- list.files(points[index])
    jsonf <- in_files[grepl("cloud_segmentation", in_files)]
    js_list <- jsonlite::read_json(sprintf("%s/%s", points[index], jsonf))
    js_list$classes <- list(
      list(
        name = "Clear",
        description = "All clear pixels, i.e. without cloud contamination or cloud shadows.",
        colour = c(255,255,255,0),
        user_colour = c(0,255,255,70)
      ),
      list(
        name = "TkCloud - center",
        description = "All cloudy pixels covered by thick clouds at borders.",
        colour = c(255, 255, 0, 70)
      ),
      list(
        name = "TkCloud - edges",
        description = "All cloudy pixels covered by thick clouds at the center.",
        colour = c(200, 200, 120, 70)
      ),
      list(
        name = "TnCloud - center",
        description = "Clouds that are semi-transparent, i.e. one can see land or sea surfaces through them. If a thin cloud lays over a thick cloud, please paint them with the <i>Thick Cloud</i> class.",
        colour = c(0, 255, 0, 70)
      ),
      list(
        name = "TnCloud - edges",
        description = "Clouds that are semi-transparent, i.e. one can see land or sea surfaces through them. If a thin cloud lays over a thick cloud, please paint them with the <i>Thick Cloud</i> class.",
        colour = c(50, 140, 50, 70)
      ),
      list(
        name = "cShadows - center",
        description = "All pixels contaminated by cloud shadows (not terrain shadows).",
        colour = c(255, 0, 0, 70)
      ),
      list(
        name = "cShadows - edges",
        description = "All pixels contaminated by cloud shadows (not terrain shadows).",
        colour = c(180, 50, 50, 70)
      ),
      list(
        name = "No data",
        description = "Reserved for no data pixels, e.g. pixels outside of the satellite's swath.",
        colour = c(160, 50, 50, 70)
      )
    )
    jsonlite::write_json(
      x = js_list,
      path = sprintf("%s/%s", points[index], jsonf),
      pretty = TRUE,
      auto_unbox = TRUE
    )
  }
}

folder <- "PUT HERE YOUR FOLDER"
scribble_creator(folder)


