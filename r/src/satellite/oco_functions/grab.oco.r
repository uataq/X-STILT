# ---------------------------------------------------------------------------- #
#' subroutine to grab OCO-2/3 data, given a spatial domain and a date
# ---------------------------------------------------------------------------- #
#' add view angles and azimuth angles, 09/24/2018 
#' update for lite v9 data, DW, 10/15/2018 
#' update for OCO-3 early data, DW, 06/28/2020 
# ---------------------------------------------------------------------------- #

# need to consider both situations, DW, 12/18/2020 
# 1. same oco orbit in a day, but the overpass time ranges from one hour to another
# 2. multiple oco orbits in a day -> different orbit and different overpass hour

grab.oco = function(oco.path, timestr, lon.lat, oco.ver = 'b9r', oco.fn = NULL) {

  library(ncdf4); library(dplyr)

  if (is.null(oco.fn)) {
    fn = list.files(oco.path, paste0('_', substr(timestr, 3, 8), '_'), 
                    full.names = T, recursive = T)
    if (length(fn) == 0) { cat('NO OCO files found...please check..\n'); return() }
  } else fn = oco.fn

  oco.dat = nc_open(fn)

  ## grabbing OCO-2 levels, lat, lon
  # level 1 to 20, for space-to-surface, level 20 is the bottom level
  oco.level = ncvar_get(oco.dat, 'levels')
  oco.lat   = ncvar_get(oco.dat, 'latitude')
  oco.lon   = ncvar_get(oco.dat, 'longitude')
  
  # OCO-3 orbit number since mission start
  oco.orbit = ncvar_get(oco.dat, 'Sounding/orbit')
  xco2.obs  = ncvar_get(oco.dat, 'xco2')
  xco2.obs[xco2.obs == -999999] = NA
  xco2.obs.uncert = ncvar_get(oco.dat, 'xco2_uncertainty')

  # get lat/lon corners, ndim of 4
  oco.lons = ncvar_get(oco.dat, 'vertex_longitude')
  oco.lats = ncvar_get(oco.dat, 'vertex_latitude')
  vertices = ncvar_get(oco.dat, 'vertices')
  dimnames(oco.lons) = list(vertices, 1 : length(oco.lat))
  dimnames(oco.lats) = list(vertices, 1 : length(oco.lat))
  oco.vert.df = full_join(reshape2::melt(oco.lons), reshape2::melt(oco.lats), 
                           by = c('Var1', 'Var2')) %>% 
                 dplyr::rename(vertices = Var1, indx = Var2, 
                               lons = value.x, lats = value.y)


  # Warn level being removed for lite v9 data, DW, 10/15/2018 
  if (grepl('7', oco.ver) | grepl('8', oco.ver)) wl = ncvar_get(oco.dat, 'warn_level')
  qf = ncvar_get(oco.dat, 'xco2_quality_flag')
  foot = ncvar_get(oco.dat, 'Sounding/footprint')
  psurf = ncvar_get(oco.dat, 'Retrieval/psurf')  # hpa
  aod.tot = ncvar_get(oco.dat, 'Retrieval/aod_total')  # Total Cloud+Aerosol Optical Depth
  aod.fine = ncvar_get(oco.dat, 'Retrieval/aod_sulfate') + 
              ncvar_get(oco.dat, 'Retrieval/aod_oc')    # AOD for sulfate + OC
  #ws.surf = ncvar_get(oco.dat, 'Retrieval/windspeed')  # m/s

  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  id   = as.character(ncvar_get(oco.dat, 'sounding_id'))
  sec  = ncvar_get(oco.dat, 'time')
  time = as.POSIXct(sec, origin = '1970-01-01 00:00:00', tz = 'UTC')
  hr   = as.numeric(substr(time, 12, 13))
  min  = as.numeric(substr(time, 15, 16))
  sec  = as.numeric(substr(time, 18, 19))

  # 0:Nadir, 1:Glint, 2:Target, 3: Transition, 4=Snapshot Area Map
  OM = ncvar_get(oco.dat, 'Sounding/operation_mode')
  LF = ncvar_get(oco.dat, 'Sounding/land_fraction') # >= 80%: land, < 20%: sea

  # operation modes:
  mode = rep(NULL, length(OM))
  mode[LF >= 80 & OM == 0] = 'Land_Nadir'
  mode[LF <= 20 & OM == 0] = 'Sea_Nadir'
  mode[LF >= 80 & OM == 1] = 'Land_Glint'
  mode[LF <= 20 & OM == 1] = 'Sea_Glint'
  mode[LF >= 80 & OM == 2] = 'Land_Target'
  mode[LF <= 20 & OM == 2] = 'Sea_Target'
  mode[LF >= 80 & OM == 3] = 'Land_Transition'
  mode[LF <= 20 & OM == 3] = 'Sea_Transition'

  if (4 %in% OM) {
    mode[LF >= 80 & OM == 4] = 'Land_SAM'
    mode[LF <= 20 & OM == 4] = 'Sea_SAM'
  }

  # grab zenith angles and azimuth angles, 09/24/2018 
  # solar zenith angle at the time of the measurement
  sza = ncvar_get(oco.dat, 'solar_zenith_angle')  # sounding_solar_zenith

  # zenith angle of the satellite at the time of the measurement
  oza = ncvar_get(oco.dat, 'sensor_zenith_angle') # sounding_zenith

  # solar azimuth angle at the time of the measurement; degrees East of North
  # same as 'sounding_solar_azimuth' in full product
  saa = ncvar_get(oco.dat, 'Sounding/solar_azimuth_angle') 

  # azimuth angle of the satellite at the time of the measurement
  # same as 'sounding_azimuth' in full product
  oaa = ncvar_get(oco.dat, 'Sounding/sensor_azimuth_angle')

  # Angular distance from viewing along the perfect glint direction
  ga = ncvar_get(oco.dat, 'Sounding/glint_angle')  # degrees

  # Airmass, computed as 1/cos(solar_zenith_angle) + 1/cos(sensor_zenith_angle)
  air.mass = ncvar_get(oco.dat, 'Sounding/airmass')  # degrees
  

  # combine all info into a data frame
  obs.all = data.frame(id = as.numeric(id), orbit = as.numeric(oco.orbit), 
                        time = as.character(time), hr, min, sec, 
                        lat = as.numeric(oco.lat), lon = as.numeric(oco.lon), 
                        foot = as.numeric(foot), qf = as.numeric(qf), 
                        psurf = as.numeric(psurf), xco2 = as.numeric(xco2.obs), 
                        xco2.uncert = as.numeric(xco2.obs.uncert),
                        mode = as.character(mode), aod.tot = as.numeric(aod.tot), 
                        aod.fine = as.numeric(aod.fine), sza, oza, saa, oaa, ga, 
                        air.mass, stringsAsFactors = F)

  if (grepl('7', oco.ver) | grepl('8', oco.ver)) 
    obs.all = cbind(obs.all, wl = as.numeric(wl))

  # add lat/lon corners
  obs.vert = obs.all %>% mutate(indx = as.numeric(rownames(obs.all))) %>% 
                         left_join(oco.vert.df, by = 'indx')
                          
  # select regions, lon.lat: c(minlon, maxlon, minlat, maxlat) with buffer of 0.01deg
  obs = obs.vert %>% filter(lat >= lon.lat$minlat - 0.01, 
                            lat <= lon.lat$maxlat + 0.01, 
                            lon >= lon.lat$minlon - 0.01, 
                            lon <= lon.lat$maxlon + 0.01) %>% 
                      mutate(timestr = substr(id, 1, 10))

  # if there are multiple orbits, further subset `obs` based on `timestr`, DW, 12/18/2020
  if ( length(unique(obs$orbit)) > 1 & nchar(timestr) == 10 ) {
    cat(paste('grab.oco(): detected multiple orbits in one day, chose soundings close to', timestr, '\n'))
    obs = obs[obs$timestr == timestr, ]
  } # end if

  sel.mode = unique(obs$mode); cat('Operational Modes:', unique(sel.mode), '\n')
  nc_close(oco.dat)
  
  return(obs)
}
# end of subroutinr

