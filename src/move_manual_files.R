#' Move manual.png to a single folder
#' 
#' Google Drive sucks for managing a large dataset so we decide to do 
#' everything locally (I had to buy other hard disks of 5 Tb T-T ). 
#' This code helps the labelers to copy only the manual labeling and avoid 
#' upload too many files.
#' @author Cesar Aybar

download_manual_labeling <- function(label_path, result_path){
  full_points <- list.files(label_path,full.names = TRUE)
  message(sprintf("Hay %s puntos", length(full_points)))
  if (!all(grepl("point_", basename(full_points)))) {
    stop(
      "Existen carpetas o archivos que no empiezan con 'point_', eliminar primero",
      " todos los archivos/carpetas que no esten relacionado con cloudSEN12."
    )
  }
  # Crear carpeta resultados
  dir.create(result_path, showWarnings = FALSE)
  
  for(index in seq_along(full_points)) {
    manual_files <- list.files(
      path = full_points[index],
      pattern = "manual\\.png",
      full.names = TRUE,
      recursive = TRUE
    )
    
    if (length(manual_files) != 5) {
      stop(
        sprintf(
          "El punto %s tiene solo %s cuando deberia tener %s",
          basename(full_points[index]),
          length(manual_files),
          "5 archivos con el nombre manual.png. REVISAR!!"
        )
      )
    }
    
    new_name <- sprintf(
      "%s/%s__%s__manual.png",
      result_path,
      basename(full_points[index]), 
      basename(dirname(dirname(manual_files)))
    )
    
    file.copy(
      manual_files,
      new_name,
      overwrite = TRUE,
    )
    message(sprintf("Punto %s copiado con exito!", basename(full_points[index])))
  }
}

label_path <- "/home/csaybar/Desktop/demo/"
result_path <- "/home/csaybar/Desktop/result/"

download_manual_labeling(label_path, result_path)
