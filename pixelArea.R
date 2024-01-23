pixelArea <- function(r, # SpatRaster
                      type = "mean", # can also be "centroid"
                      unit = "m",  # can also be "km"
                      mask = TRUE,  # to consider only areas of non-NA pixels
                      map = TRUE) {
  # version 1.0 (23 Jan 2024)
  # https://modtools.wordpress.com/2024/01/23/actual-pixel-sizes-of-unprojected-raster-maps/
  
  stopifnot(inherits(r, "SpatRaster"),
            type %in% c("mean", "centroid"))
  
  r_size <- terra::cellSize(r, unit = unit, mask = mask)
  if (map) terra::plot(r_size, main = paste0("Pixel area (", unit, "2)"))
  areas <- terra::values(r_size, mat = FALSE, dataframe = FALSE, na.rm = FALSE)  # na.rm must be FALSE for areas[centr_pix] to be correct
  
  if (type == "mean") {
    cat(paste0("Mean pixel area (", unit, "2):\n"))
    return(mean(areas, na.rm = TRUE))
  }
  
  if (type == "centroid") {
    r_pol <- terra::as.polygons(r * 0, aggregate = TRUE)
    centr <- terra::centroids(r_pol)
    if (map) {
      if (!mask) terra::plot(r_pol, lwd = 0.2, add = TRUE)
      terra::plot(centr, pch = 4, col = "blue", add = TRUE)
    }
    
    centr_pix <- terra::cellFromXY(r, terra::crds(centr))
    cat(paste0("Centroid pixel area (", unit, "2):\n"))
    return(areas[centr_pix])
  }
}
