projectCentered <- function(x,  # SpatRaster or SpatVector
                            proj = "laea"  # +proj parameter for PROJ.4 string
                            ) {
  # version 1.0, by A. Marcia Barbosa
  stopifnot(inherits(x, "SpatRaster") || inherits(x, "SpatVector"))
  if (terra::crs(x) == "" || is.na(terra::crs(x))) {
    stop("'x' does not have a CRS to project")
  }
  if (inherits(x, "SpatRaster")) {
    a <- terra::as.polygons(x * 0)
  }
  if (inherits(x, "SpatVector") && length(x) > 1) {
    a <- terra::aggregate(x)
  } else {
    a <- x
  }
  if (!isTRUE(terra::is.lonlat(a))) {
    a <- terra::project(a, "EPSG:4326")
  }
  centr <- terra::crds(terra::centroids(a))
  proj <- tolower(proj)
  proj <- sub("+proj=", "", proj, fixed = TRUE)  # in case user input prefix
  CRS <- paste(
    paste0("+proj=", proj),
    paste0("+lat_0=", centr[2]),
    paste0("+lon_0=", centr[1])
  )
  tryCatch(terra::project(x, CRS),
           error = function(e)
             stop("'proj' is invalid"))
}
