# subroutine to grab footprint from STILT output and convert to data.frame
# DW, 07/10/2018
# remove footpath,

grab.foot <- function(stilt.ver, footfile, foot.sig = 1E-6, lon.lat = NULL){

  if (stilt.ver == 1) {
    library(ncdf4); library(reshape)

    footdat <- nc_open(footfile)
    foot <- ncvar_get(footdat, 'footprint')
    lat  <- ncvar_get(footdat, 'Lat') + 1/240  # move to centered lat lon
    lon  <- ncvar_get(footdat, 'Lon') + 1/240
    dimnames(foot) <- list(lat, lon)

    melt.foot <- reshape2::melt(foot)
    colnames(melt.foot) <- list('lat', 'lon', 'foot')
    nc_close(footdat)

  } else {
    library(raster)

    # footprint generated from Ben's codes
    foot <- raster(footfile)
    melt.foot <- raster::as.data.frame(x = foot, xy = T)
    colnames(melt.foot) <- list('lon', 'lat', 'foot')
  }

  # select footprint domain and footprint values
  sel.foot <- melt.foot %>% filter(foot >= foot.sig)

  if (!is.null(lon.lat)) {
    sel.foot <- sel.foot %>% filter(
      lon >= lon.lat$minlon & lon <= lon.lat$maxlon &
      lat >= lon.lat$minlat & lat <= lon.lat$maxlat)
  }

  return(sel.foot)
}
# end of subroutine
