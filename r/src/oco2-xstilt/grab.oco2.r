#' subroutine to grab OCO-2 data, given a spatial domain and a date
#' @author: Dien Wu

#' @updates: 
#' add view angles and azimuth angles, 09/24/2018 
#' update for lite v9 data, DW, 10/15/2018 

grab.oco2 <- function(ocopath, timestr, lon.lat, oco2.ver = 'b9r'){

  library(ncdf4)
  ocofile <- list.files(pattern = substr(timestr, 3, 8), path = ocopath)
  if (length(ocofile) == 0) {
    cat('NO OCO2 files found...please check..\n'); return()
  }
  ocodat <- nc_open(file.path(ocopath, ocofile))

  ## grabbing OCO-2 levels, lat, lon
  # level 1 to 20, for space-to-surface, level 20 is the bottom level
  oco.level <- ncvar_get(ocodat, 'levels')
  oco.lat   <- ncvar_get(ocodat, 'latitude')
  oco.lon   <- ncvar_get(ocodat, 'longitude')
  xco2.obs  <- ncvar_get(ocodat, 'xco2')
  xco2.obs[xco2.obs == -999999] <- NA
  xco2.obs.uncert <- ncvar_get(ocodat, 'xco2_uncertainty')

  # Warn level being removed for lite v9 data, DW, 10/15/2018 
  if (oco2.ver != 'b9r') wl <- ncvar_get(ocodat, 'warn_level')
  qf    <- ncvar_get(ocodat, 'xco2_quality_flag')
  foot  <- ncvar_get(ocodat, 'Sounding/footprint')
  psurf <- ncvar_get(ocodat, 'Retrieval/psurf')  # hpa
  #ws.surf <- ncvar_get(ocodat, 'Retrieval/windspeed')  # m/s

  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  id   <- as.character(ncvar_get(ocodat, 'sounding_id'))
  sec  <- ncvar_get(ocodat, 'time')
  time <- as.POSIXct(sec, origin = '1970-01-01 00:00:00', tz = 'UTC')
  hr   <- as.numeric(substr(time, 12, 13))
  min  <- as.numeric(substr(time, 15, 16))
  sec  <- as.numeric(substr(time, 18, 19))

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

  # grab zenith angles and azimuth angles, 09/24/2018 
  # solar zenith angle at the time of the measurement
  sza <- ncvar_get(ocodat, 'solar_zenith_angle')  # sounding_solar_zenith

  # zenith angle of the satellite at the time of the measurement
  oza <- ncvar_get(ocodat, 'sensor_zenith_angle') # sounding_zenith

  # solar azimuth angle at the time of the measurement; degrees East of North
  # same as 'sounding_solar_azimuth' in full product
  saa <- ncvar_get(ocodat, 'Sounding/solar_azimuth_angle') 

  # azimuth angle of the satellite at the time of the measurement
  # same as 'sounding_azimuth' in full product
  oaa <- ncvar_get(ocodat, 'Sounding/sensor_azimuth_angle')

  # Angular distance from viewing along the perfect glint direction
  ga <- ncvar_get(ocodat, 'Sounding/glint_angle')  # degrees

  # Airmass, computed as 1/cos(solar_zenith_angle) + 1/cos(sensor_zenith_angle)
  air.mass <- ncvar_get(ocodat, 'Sounding/airmass')  # degrees
  
  obs.all <- data.frame(id = as.numeric(id), time = as.character(time), 
                        hr, min, sec, lat = as.numeric(oco.lat), 
                        lon = as.numeric(oco.lon), foot = as.numeric(foot), 
                        qf = as.numeric(qf), 
                        xco2 = as.numeric(xco2.obs), 
                        xco2.uncert = as.numeric(xco2.obs.uncert),
                        mode = as.character(mode), sza, oza, saa, oaa, ga, 
                        air.mass, stringsAsFactors = F)
  if (oco2.ver != 'b9r') obs.all <- cbind(obs.all, wl = as.numeric(wl))

  # select regions, lon.lat: c(minlon, maxlon, minlat, maxlat)
  obs <- obs.all %>% 
    filter(lat >= lon.lat$minlat & lat <= lon.lat$maxlat &
           lon >= lon.lat$minlon & lon <= lon.lat$maxlon) %>% 
    mutate(timestr = substr(id, 1, 10))

  sel.mode <- unique(obs$mode)
  #cat('Operational Modes:', unique(sel.mode), '\n\n')

  nc_close(ocodat)
  return(obs)
}
# end of subroutinr