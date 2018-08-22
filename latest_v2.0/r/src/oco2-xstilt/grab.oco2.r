# subroutine to grab OCO-2 data, given a spatial domain and a date
# DW

grab.oco2 <- function(ocopath, timestr, lon.lat){

  library(ncdf4)
  ocofile <- list.files(pattern = substr(timestr, 3, 8), path = ocopath)
  ocodat <- nc_open(file.path(ocopath, ocofile))

  ## grabbing OCO-2 levels, lat, lon
  # level 1 to 20, for space-to-surface, level 20 is the bottom level
  oco.level <- ncvar_get(ocodat, "levels")
  oco.lat   <- ncvar_get(ocodat, "latitude")
  oco.lon   <- ncvar_get(ocodat, "longitude")
  xco2.obs  <- ncvar_get(ocodat, "xco2")
  xco2.obs[xco2.obs == -999999] <- NA
  xco2.obs.uncert <- ncvar_get(ocodat, "xco2_uncertainty")

  wl <- ncvar_get(ocodat, "warn_level")
  qf <- ncvar_get(ocodat, "xco2_quality_flag")
  foot <- ncvar_get(ocodat, "Sounding/footprint")

  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  id <- as.character(ncvar_get(ocodat, "sounding_id"))
  sec <- ncvar_get(ocodat, "time")
  time <- as.POSIXct(sec, origin = "1970-01-01 00:00:00", tz = "UTC")

  # 0:Nadir, 1:Glint, 2:Target, 3: Transition
  OM <- ncvar_get(ocodat, 'Sounding/operation_mode')
  LF <- ncvar_get(ocodat, 'Sounding/land_fraction') # > 80%: land, < 20%: sea

  # operation modes:
  mode <- rep(NULL, length(OM))
  mode[LF > 80 & OM == 0] <- 'Land_Nadir'
  mode[LF < 20 & OM == 0] <- 'Sea_Nadir'
  mode[LF > 80 & OM == 1] <- 'Land_Glint'
  mode[LF < 20 & OM == 1] <- 'Sea_Glint'
  mode[LF > 80 & OM == 2] <- 'Land_Target'
  mode[LF < 20 & OM == 2] <- 'Sea_Target'
  mode[LF > 80 & OM == 3] <- 'Land_Transition'
  mode[LF < 20 & OM == 3] <- 'Sea_Transition'

  # select regions, lon.lat: c(minlon, maxlon, minlat, maxlat)
  region.index <- oco.lat >= lon.lat$minlat & oco.lat <= lon.lat$maxlat &
                  oco.lon >= lon.lat$minlon & oco.lon <= lon.lat$maxlon

  # get OCO-2 overpassing time
  sel.id   <- as.numeric(id[region.index])
  sel.time <- as.character(time[region.index])
  sel.hr   <- as.numeric(substr(sel.time, 12, 13))
  sel.min  <- as.numeric(substr(sel.time, 15, 16))
  sel.sec  <- as.numeric(substr(sel.time, 18, 19))
  sel.foot <- as.numeric(foot[region.index])
  sel.wl   <- as.numeric(wl[region.index])
  sel.qf   <- as.numeric(qf[region.index])
  sel.lat  <- as.numeric(oco.lat[region.index])
  sel.lon  <- as.numeric(oco.lon[region.index])
  sel.xco2.obs <- as.numeric(xco2.obs[region.index])
  sel.xco2.obs.uncert <- as.numeric(xco2.obs.uncert[region.index])
  sel.mode <- mode[region.index]
  cat('Operational Modes:', unique(sel.mode), '\n\n')

  obs.all <- data.frame(id = sel.id, time = sel.time, hr = sel.hr,
    min = sel.min, sec = sel.sec, lat = sel.lat, lon = sel.lon, foot = sel.foot,
    xco2 = sel.xco2.obs, xco2.uncert = sel.xco2.obs.uncert, wl = sel.wl,
    qf = sel.qf)

  return(obs.all)
}
