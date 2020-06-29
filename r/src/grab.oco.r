#' subroutine to grab OCO-2/3 data, given a spatial domain and a date
#' @author: Dien Wu

#' @updates: 
#' add view angles and azimuth angles, 09/24/2018 
#' update for lite v9 data, DW, 10/15/2018 
#' update for OCO-3 early data, DW, 06/28/2020 

grab.oco <- function(oco.path, timestr, lon.lat, oco.ver = 'b9r'){

  library(ncdf4)
  oco.file <- list.files(pattern = substr(timestr, 3, 8), path = oco.path)
  if (length(oco.file) == 0) {
    cat('NO OCO files found...please check..\n'); return()
  }
  oco.dat <- nc_open(file.path(oco.path, oco.file))

  ## grabbing OCO-2 levels, lat, lon
  # level 1 to 20, for space-to-surface, level 20 is the bottom level
  oco.level <- ncvar_get(oco.dat, 'levels')
  oco.lat   <- ncvar_get(oco.dat, 'latitude')
  oco.lon   <- ncvar_get(oco.dat, 'longitude')
  xco2.obs  <- ncvar_get(oco.dat, 'xco2')
  xco2.obs[xco2.obs == -999999] <- NA
  xco2.obs.uncert <- ncvar_get(oco.dat, 'xco2_uncertainty')

  # Warn level being removed for lite v9 data, DW, 10/15/2018 
  if (grepl('7', oco.ver) | grepl('8', oco.ver)) wl <- ncvar_get(oco.dat, 'warn_level')
  qf    <- ncvar_get(oco.dat, 'xco2_quality_flag')
  foot  <- ncvar_get(oco.dat, 'Sounding/footprint')
  psurf <- ncvar_get(oco.dat, 'Retrieval/psurf')  # hpa
  #ws.surf <- ncvar_get(oco.dat, 'Retrieval/windspeed')  # m/s

  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  id   <- as.character(ncvar_get(oco.dat, 'sounding_id'))
  sec  <- ncvar_get(oco.dat, 'time')
  time <- as.POSIXct(sec, origin = '1970-01-01 00:00:00', tz = 'UTC')
  hr   <- as.numeric(substr(time, 12, 13))
  min  <- as.numeric(substr(time, 15, 16))
  sec  <- as.numeric(substr(time, 18, 19))

  # 0:Nadir, 1:Glint, 2:Target, 3: Transition
  OM <- ncvar_get(oco.dat, 'Sounding/operation_mode')
  LF <- ncvar_get(oco.dat, 'Sounding/land_fraction') # > 80%: land, < 20%: sea

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

  # grab zenith angles and azimuth angles, 09/24/2018 
  # solar zenith angle at the time of the measurement
  sza <- ncvar_get(oco.dat, 'solar_zenith_angle')  # sounding_solar_zenith

  # zenith angle of the satellite at the time of the measurement
  oza <- ncvar_get(oco.dat, 'sensor_zenith_angle') # sounding_zenith

  # solar azimuth angle at the time of the measurement; degrees East of North
  # same as 'sounding_solar_azimuth' in full product
  saa <- ncvar_get(oco.dat, 'Sounding/solar_azimuth_angle') 

  # azimuth angle of the satellite at the time of the measurement
  # same as 'sounding_azimuth' in full product
  oaa <- ncvar_get(oco.dat, 'Sounding/sensor_azimuth_angle')

  # Angular distance from viewing along the perfect glint direction
  ga <- ncvar_get(oco.dat, 'Sounding/glint_angle')  # degrees

  # Airmass, computed as 1/cos(solar_zenith_angle) + 1/cos(sensor_zenith_angle)
  air.mass <- ncvar_get(oco.dat, 'Sounding/airmass')  # degrees
  
  obs.all <- data.frame(id = as.numeric(id), time = as.character(time), 
                        hr, min, sec, lat = as.numeric(oco.lat), 
                        lon = as.numeric(oco.lon), foot = as.numeric(foot), 
                        qf = as.numeric(qf), 
                        xco2 = as.numeric(xco2.obs), 
                        xco2.uncert = as.numeric(xco2.obs.uncert),
                        mode = as.character(mode), sza, oza, saa, oaa, ga, 
                        air.mass, stringsAsFactors = F)
  if (grepl('7', oco.ver) | grepl('8', oco.ver)) 
      obs.all <- cbind(obs.all, wl = as.numeric(wl))

  # select regions, lon.lat: c(minlon, maxlon, minlat, maxlat)
  obs <- obs.all %>% filter(lat >= lon.lat$minlat & lat <= lon.lat$maxlat &
                            lon >= lon.lat$minlon & lon <= lon.lat$maxlon) %>% 
                     mutate(timestr = substr(id, 1, 10))

  sel.mode <- unique(obs$mode)
  #cat('Operational Modes:', unique(sel.mode), '\n\n')

  nc_close(oco.dat)
  return(obs)
}
# end of subroutinr