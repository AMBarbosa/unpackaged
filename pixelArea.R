pixelArea <- function(rast, # SpatRaster
                      type = "mean", # can also be "centroid"
                      unit = "m",  # can also be "km"
                      mask = TRUE,  # to use only non-NA pixels
                      map = TRUE) {

  # by A. Marcia Barbosa (https://modtools.wordpress.com/)
  # version 1.4 (17 May 2024)

  r <- rast

  stopifnot(inherits(r, "SpatRaster"),
            type %in% c("mean", "centroid"))

  r_size <- terra::cellSize(r, unit = unit, mask = mask)
  if (map) terra::plot(r_size, main = paste0("Pixel area (", unit, "2)"))
  areas <- terra::values(r_size, mat = FALSE, dataframe = FALSE, na.rm = FALSE)  # na.rm must be FALSE for areas[centr_pix] to be correct

  if (type == "mean") {
    out <- mean(areas, na.rm = TRUE)
    cat(paste0("Mean pixel area (", unit, "2):\n", out, "\n\n"))
    return(out)
  }

  if (type == "centroid") {
    r_pol <- terra::as.polygons(r * 0, aggregate = TRUE)
    centr <- terra::centroids(r_pol)
    if (map) {
      if (!mask) terra::plot(r_pol, lwd = 0.2, add = TRUE)
      terra::plot(centr, pch = 4, col = "blue", add = TRUE)
    }

    centr_pix <- terra::cellFromXY(r, terra::crds(centr))
    out <- areas[centr_pix]
    if (!is.finite(out)) message("The centroid of your region may not have a pixel value; consider using mask=FALSE, or type = 'mean'.")
    cat(paste0("Centroid pixel area (", unit, "2):\n", out, "\n\n"))
    return(out)
  }
}
