## functions related to extracting OCO-2/3 data
# including find.overpass(), get.site.track(), grab.oco(), grab.sif(), get.oco.info()



# ---------------------------------------------------------------------------- #
# function that search for any OCO-2/3 overpases for a given region, DW, 05/15/2017
#' @param date.range date range, c(YYYYMMDD_1, YYYYMMDD_2)
#' @param oco.path default in lin-group
#' @param oco.ver needs to contain version number
#' @param lon.lat a data frame containining c(minlat, maxlat, minlon, maxlon), 
#'                please follow THIS ORDER!!!
#' @param urbanTF for searching soundings near urban region, DW, 06/15/2018
# ---------------------------------------------------------------------------- #
# add count for good quality data, DW, 12/20/2017
# update for v9 data, DW, 10/19/2018 
# drop the scientific notation for sounding ID, DR, DW, 09/04/2019
# update for OCO-3 Vearly data, DW, 06/29/2020 
# remove warn levels, DW, 07/01/2020 
# ---------------------------------------------------------------------------- #

find.overpass <- function(date.range, lon.lat, oco.ver = 'b9r', oco.path, 
                          urbanTF = F, dlon = 0.5, dlat = 0.5){ 

    library(geosphere); library(ncdf4); library(dplyr)

    # path and filename for storing OCO-2 info
    all.file <- list.files(pattern = 'LtCO2_', path = oco.path)

    # get rid of some characters
    file.info <- gsub('.nc4', '', all.file)
    file.info <- strsplit.to.df(file.info)
    all.timestr <- file.info$V3

    if (grepl('7', oco.ver)) {   # for version 7r
      all.timestr[nchar(all.timestr) == 6] <-
        paste0('20', all.timestr[nchar(all.timestr) == 6])
      oco.file <- all.file[all.timestr >= date.range[1] & all.timestr <= date.range[2]]
      timestr <- all.timestr[all.timestr >= date.range[1] & all.timestr <= date.range[2]]

    } else if (!grepl('7', oco.ver)) {   # for version 8, 9r or even for OCO-3
      SEL.day <- all.timestr >= substr(date.range[1], 3, 8) &
                 all.timestr <= substr(date.range[2], 3, 8)
      oco.file <- all.file[SEL.day]
      timestr <- paste0('20', all.timestr[SEL.day])
    } 

    # loop over each overpass
    result <- NULL
    for (f in 1 : length(oco.file)) {
      
      if (f %% 25 == 0)
      cat(paste('#-- ', signif(f / length(oco.file), 3) * 100, '% SEARCHED --#\n'))
      
      dat <- nc_open(file.path(oco.path, oco.file[f]))

      ## grabbing OCO-2 levels, lat, lon
      oco.lev <- ncvar_get(dat, 'levels')
      oco.lat <- ncvar_get(dat, 'latitude')
      oco.lon <- ncvar_get(dat, 'longitude')
      xco2 <- ncvar_get(dat, 'xco2'); xco2[xco2 == -999999] <- NA
      qf <- ncvar_get(dat, 'xco2_quality_flag')

      # drop the scientific notation for sounding ID, DR, DW, 09/04/2019
      id <- format(ncvar_get(dat, 'sounding_id'), scientific = F)

      # 0=Nadir, 1=Glint, 2=Target, 3=Transition, 4=Snapshot Area Map
      mode <- ncvar_get(dat, 'Sounding/operation_mode')

      obs <- data.frame(lat = as.numeric(oco.lat), lon = as.numeric(oco.lon), 
                        qf = as.numeric(qf), xco2 = as.numeric(xco2), 
                        timestr = as.numeric(substr(id, 1, 10)), mode = mode,
                        stringsAsFactors = F) 
      
      obs <- obs %>% filter(lat >= lon.lat$minlat & lat <= lon.lat$maxlat &
                            lon >= lon.lat$minlon & lon <= lon.lat$maxlon) 

      
      # count overpasses
      tot.count <- as.numeric(nrow(obs))
      qf.count  <- as.numeric(nrow(obs %>% filter(qf == 0)))
      sam.count <- as.numeric(nrow(obs %>% filter(mode == 4)))
      nc_close(dat)
        
      # store results if there are soundings over
      if (tot.count > 0) {
        uni.timestr <- unique(obs$timestr)

        tmp <- data.frame(timestr = uni.timestr, tot.count = tot.count,
                          sam.count = sam.count, qf.count = qf.count, 
                          stringsAsFactors = F)

        # also search for soundings near city center,
        if (urbanTF) {
          urban.dat <- obs %>% filter(lon >= (lon.lat$citylon - dlon),
                                      lon <= (lon.lat$citylon + dlon),
                                      lat >= (lon.lat$citylat - dlat),
                                      lat <= (lon.lat$citylat + dlat))
          tot.urban.count <- as.numeric(nrow(urban.dat))
          qf.urban.count  <- as.numeric(nrow(urban.dat %>% filter(qf == 0)))
          sam.urban.count <- as.numeric(nrow(urban.dat %>% filter(qf == 0, mode == 4)))

          # combine 
          tmp <- tmp %>% mutate(tot.urban.count = tot.urban.count, 
                                sam.urban.count = sam.urban.count,
                                qf.urban.count = qf.urban.count)
        } # end if urbanTF

        result <- rbind(result, tmp)
      } else next  # if no data, jump to next file
      # end if tot.count

    }  # end for f

    result <- result[order(result$timestr),]
    return(result)
}


# ---------------------------------------------------------------------------- #
#' subroutine to grab OCO-2/3 data, given a spatial domain and a date
# ---------------------------------------------------------------------------- #
#' add view angles and azimuth angles, 09/24/2018 
#' update for lite v9 data, DW, 10/15/2018 
#' update for OCO-3 early data, DW, 06/28/2020 
# ---------------------------------------------------------------------------- #

grab.oco <- function(oco.path, timestr, lon.lat, oco.ver = 'b9r'){

  library(ncdf4)
  oco.file <- list.files(oco.path, paste0('_', substr(timestr, 3, 8), '_'))
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

  # get lat/lon corners, ndim of 4
  oco.lons  <- ncvar_get(oco.dat, 'vertex_longitude')
  oco.lats  <- ncvar_get(oco.dat, 'vertex_latitude')
  vertices  <- ncvar_get(oco.dat, 'vertices')
  dimnames(oco.lons) <- list(vertices, 1 : length(oco.lat))
  dimnames(oco.lats) <- list(vertices, 1 : length(oco.lat))
  oco.vert.df <- full_join(melt(oco.lons), melt(oco.lats), by = c('Var1', 'Var2')) %>% 
                 rename(vertices = Var1, indx = Var2, lons = value.x, lats = value.y)


  # Warn level being removed for lite v9 data, DW, 10/15/2018 
  if (grepl('7', oco.ver) | grepl('8', oco.ver)) wl <- ncvar_get(oco.dat, 'warn_level')
  qf    <- ncvar_get(oco.dat, 'xco2_quality_flag')
  foot  <- ncvar_get(oco.dat, 'Sounding/footprint')
  psurf <- ncvar_get(oco.dat, 'Retrieval/psurf')  # hpa
  aod.tot <- ncvar_get(oco.dat, 'Retrieval/aod_total')  # Total Cloud+Aerosol Optical Depth
  aod.fine <- ncvar_get(oco.dat, 'Retrieval/aod_sulfate') + 
              ncvar_get(oco.dat, 'Retrieval/aod_oc')    # AOD for sulfate + OC
  #ws.surf <- ncvar_get(oco.dat, 'Retrieval/windspeed')  # m/s

  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  id   <- as.character(ncvar_get(oco.dat, 'sounding_id'))
  sec  <- ncvar_get(oco.dat, 'time')
  time <- as.POSIXct(sec, origin = '1970-01-01 00:00:00', tz = 'UTC')
  hr   <- as.numeric(substr(time, 12, 13))
  min  <- as.numeric(substr(time, 15, 16))
  sec  <- as.numeric(substr(time, 18, 19))

  # 0:Nadir, 1:Glint, 2:Target, 3: Transition, 4=Snapshot Area Map
  OM <- ncvar_get(oco.dat, 'Sounding/operation_mode')
  LF <- ncvar_get(oco.dat, 'Sounding/land_fraction') # >= 80%: land, < 20%: sea

  # operation modes:
  mode <- rep(NULL, length(OM))
  mode[LF >= 80 & OM == 0] <- 'Land_Nadir'
  mode[LF <= 20 & OM == 0] <- 'Sea_Nadir'
  mode[LF >= 80 & OM == 1] <- 'Land_Glint'
  mode[LF <= 20 & OM == 1] <- 'Sea_Glint'
  mode[LF >= 80 & OM == 2] <- 'Land_Target'
  mode[LF <= 20 & OM == 2] <- 'Sea_Target'
  mode[LF >= 80 & OM == 3] <- 'Land_Transition'
  mode[LF <= 20 & OM == 3] <- 'Sea_Transition'

  if (4 %in% OM) {
    mode[LF >= 80 & OM == 4] <- 'Land_SAM'
    mode[LF <= 20 & OM == 4] <- 'Sea_SAM'
  }

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
  

  # combine all info into a data frame
  obs.all <- data.frame(id = as.numeric(id), time = as.character(time), 
                        hr, min, sec, lat = as.numeric(oco.lat), 
                        lon = as.numeric(oco.lon), foot = as.numeric(foot), 
                        qf = as.numeric(qf), psurf = as.numeric(psurf), 
                        xco2 = as.numeric(xco2.obs), 
                        xco2.uncert = as.numeric(xco2.obs.uncert),
                        mode = as.character(mode), aod.tot = as.numeric(aod.tot), 
                        aod.fine = as.numeric(aod.fine), sza, oza, saa, oaa, ga, 
                        air.mass, stringsAsFactors = F)

  if (grepl('7', oco.ver) | grepl('8', oco.ver)) 
    obs.all <- cbind(obs.all, wl = as.numeric(wl))

  # add lat/lon corners
  obs.vert <- obs.all %>% mutate(indx = as.numeric(rownames(obs.all))) %>% 
                          left_join(oco.vert.df, by = 'indx')
                          
  # select regions, lon.lat: c(minlon, maxlon, minlat, maxlat) with buffer of 0.01deg
  obs <- obs.vert %>% filter(lat >= lon.lat$minlat - 0.01, 
                             lat <= lon.lat$maxlat + 0.01, 
                             lon >= lon.lat$minlon - 0.01, 
                             lon <= lon.lat$maxlon + 0.01) %>% 
                      mutate(timestr = substr(id, 1, 10))

  sel.mode <- unique(obs$mode); cat('Operational Modes:', unique(sel.mode), '\n\n')

  nc_close(oco.dat)
  
  return(obs)
}
# end of subroutinr





# ---------------------------------------------------------------------------- #
# subroutine to read and select SIF data, DW 08/07/2018
# ---------------------------------------------------------------------------- #
# sif.ver same as oco.ver, e.g., V7r, V8r, V10r, VEarlyr
# ---------------------------------------------------------------------------- #
grab.sif <- function(sif.path, timestr, lon.lat, sif.ver) {

    # get Sif file name
    sif.file <- list.files(pattern = paste0('LtSIF_', substr(timestr, 3, 8)),
                           path = sif.path, full.names = T)
    oco.sensor <- substr(basename(sif.file), 1, 4)

    if (length(sif.file) == 0) {
      warnings('NO SIF file found for this OCO-2 overpass...'); return()

    } else {  # if file found...

      library(ncdf4)
      sif.dat <- nc_open(sif.file)
      
      if (oco.sensor == 'oco2' & sif.ver != 'V10r') {
        lat     <- ncvar_get(sif.dat, 'latitude')
        lon     <- ncvar_get(sif.dat, 'longitude')
        sif757  <- ncvar_get(sif.dat, 'SIF_757nm') # unit in W/m2/sr/µm
        sif771  <- ncvar_get(sif.dat, 'SIF_771nm')
        igbp    <- ncvar_get(sif.dat, 'IGBP_index')

        sif <- data.frame(timestr = as.numeric(timestr), lat = as.numeric(lat),
                          lon = as.numeric(lon), sif757 = as.numeric(sif757),
                          sif771 = as.numeric(sif771), igbp = as.numeric(igbp))

      } else {
        
        #"seconds since 1993-01-01 00:00:00 UTC"
        sec     <- ncvar_get(sif.dat, 'Geolocation/time_tai93')
        lat     <- ncvar_get(sif.dat, 'Latitude')
        lon     <- ncvar_get(sif.dat, 'Longitude')
        sif757  <- ncvar_get(sif.dat, 'Science/SIF_757nm') # unit in W/m2/sr/µm
        sif771  <- ncvar_get(sif.dat, 'Science/SIF_771nm')
        igbp    <- ncvar_get(sif.dat, 'Science/IGBP_index')

        sif <- data.frame(timestr = as.numeric(timestr), lat = as.numeric(lat),
                          lon = as.numeric(lon), sif757 = as.numeric(sif757),
                          sif771 = as.numeric(sif771), igbp = as.numeric(igbp))
      } 
     
      # select SIF in given region and scale SIF_771nm with sacling factor of
      # 1.35 (used in Luus et al., 2017) to calculate an averaged SIF
      sel.sif <- sif %>% filter(lon >= lon.lat$minlon & lon <= lon.lat$maxlon &
                                lat >= lon.lat$minlat & lat <= lon.lat$maxlat) %>%
                         mutate(avg.sif = (sif757 + sif771 * 1.35)/2)

      # assign months and seasons, only for NH for now, DW
      sel.sif <- sel.sif %>% 
                 mutate(mon = substr(timestr, 5, 6),
                        season = ifelse(mon %in% c('12', '01', '02'), 'WINTER',
                                 ifelse(mon %in% c('03', '04', '05'), 'SPRING',
                                 ifelse(mon %in% c('06', '07', '08'), 'SUMMER', 
                                                                      'FALL'))))
      nc_close(sif.dat)
      return(sel.sif)
    } # end if

} # end if subroutine




# ---------------------------------------------------------------------------- #
# function to find out the AK*pwf, apriori profiles are,
# given the lat lon time from model receptors, Dien Wu, 01/11/2017
# ---------------------------------------------------------------------------- #
### Input variables--
#' @param oco.path oco path and file for searching satellite soundings;
#' @param timestr for finding correct oco.file
#' @param recp.lat @param recp.lon numeric numbers for receptor lat/lon, one at a time
#' @param diff.td allowable thredshold for difference in lat/lon between
#'                given receptors and all satellite soundings
# ---------------------------------------------------------------------------- #
### Updates--
# replace 'timestr', 'recp.lat', 'recp.lon' with 'receptor' from 'output', DW
# minor update for using OCO-3 data, i.e., change variable names, DW, 06/28/2020
# ---------------------------------------------------------------------------- #

get.oco.info <- function(oco.path, receptor, diff.td = 1E-4){

  # grabbing OCO-2 info
  timestr  <- strftime(receptor$run_time, tz = 'UTC', format = '%Y%m%d%H')
  oco.file <- list.files(oco.path, pattern = paste0('_', substr(timestr, 3, 8), '_'))
  if (length(oco.file) == 0) stop('No OCO file found for this timestr\n')
  oco.dat  <- nc_open(file.path(oco.path, oco.file))

  # grabbing OCO-2 levels, lat, lon
  # level 1 to 20, for space-to-surface, level 20 is the bottom level
  # may need to reverse later
  oco.lev <- ncvar_get(oco.dat, 'levels')
  oco.lat <- ncvar_get(oco.dat, 'latitude')
  oco.lon <- ncvar_get(oco.dat, 'longitude')
  
  # grabbing sounding ID for STILT receptors
  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  oco.id  <- as.character(ncvar_get(oco.dat, 'sounding_id'))

  # locate the OCO2 data using lat, lon, when diff are both the smallest
  diff.lat <- abs(oco.lat - receptor$lati)
  diff.lon <- abs(oco.lon - receptor$long)

  # try to find the closest sounding lat/lon,
  # given receptor lat/lon and allowable difference thredshold, 'diff.td'
  # choose only if lat/lon indices are the same
  loc.index <- intersect(which(diff.lat < diff.td), which(diff.lon < diff.td))

  # cannot find the sounding according to receptor lat/lon
  # if so, loose 'diff.td', or check OCO-2 version, or input lat/lon
  if(length(loc.index) != 1){

    cat('get.oco.info(): cannot find the receptor lat/lon from OCO file...')
    return()

  } else {

    # if one OCO-2 sounding found for a given receptor lat/lon --
    # return the oco lat, lon, ak, pwf, apriori, profiles
    find.lat <- oco.lat[loc.index]
    find.lon <- oco.lon[loc.index]
    find.id  <- oco.id [loc.index]

    ## grab column co2, averaging kernel, pressure weight and prior CO2 profiles
    ## dimensions--[levels, soundingID]
    ap <- ncvar_get(oco.dat, 'co2_profile_apriori')[, loc.index] # in ppm
    pwf <- ncvar_get(oco.dat, 'pressure_weight')[, loc.index]  # pwf
    pres <- ncvar_get(oco.dat, 'pressure_levels')[, loc.index] # press in hPa

    # normalized averaging kernel (unitless)
    ak.norm <- ncvar_get(oco.dat, 'xco2_averaging_kernel')[, loc.index]

    ## dimensions--[soundingID]
    xco2 <- ncvar_get(oco.dat, 'xco2')[loc.index]
    grdhgt <- ncvar_get(oco.dat, 'Sounding/altitude')[loc.index] # mASL
    xco2.uncert <- ncvar_get(oco.dat, 'xco2_uncertainty')[loc.index] # ret err

    # total column water vapor in kg m-2, convert to mol m-2
    xh2o <- ncvar_get(oco.dat, 'Retrieval/tcwv')[loc.index] * 1E3 / 18  

    # satellite footprint, 1-8
    footprint <- ncvar_get(oco.dat, 'Sounding/footprint')[loc.index]
    psfc  <- ncvar_get(oco.dat, 'Retrieval/psurf')[loc.index] # sfc pressure

    # check whether is missing data
    ap[ap == -999999] <- NA
    pwf[pwf == -999999] <- NA
    xco2[xco2 == -999999] <- NA
    pres[pres == -999999] <- NA
    psfc[psfc == -999999] <- NA
    grdhgt[grdhgt == -999999] <- NA
    ak.norm[ak.norm == -999999] <- NA
    footprint[footprint == -999999] <- NA
    xco2.uncert[xco2.uncert == -999999] <- NA

    # assign vertical dimnames
    attributes(ap)$names      <- oco.lev
    attributes(pwf)$names     <- oco.lev
    attributes(pres)$names    <- oco.lev
    attributes(ak.norm)$names <- oco.lev

    ### combine all OCO-2 vertical profiles and other 1D variables
    all.info <- list(oco.id = find.id, oco.lat = find.lat, oco.lon = find.lon, 
                     ak.norm = ak.norm, pwf = pwf, pres = pres, ap = ap, 
                     oco.grdhgt = grdhgt, oco.psfc = psfc, oco.foot = footprint,
                     oco.xco2 = xco2, oco.xco2.uncert = xco2.uncert, oco.xh2o = xh2o)
    
    nc_close(oco.dat)
    all.info      # return both profiles and other retrivals
  }
  
} # end of subroutine
