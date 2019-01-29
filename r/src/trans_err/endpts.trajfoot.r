# subroutine to readin CarbonTracker-NeatRealTime for OCO-2 project
# and derive background CO2 for each trajec (with/without transport errors)
# by Dien Wu, 04/12/2017

# fix 1, 04/19/2017, DW, add AK & PW weighting to background CO2
# fix 2, remove fix1, 04/20/2017, DW, not to weight background CO2 for now
# generalizt to merge with Ben's code, DW, 07/23/2018
# read bio fluxes as raster layer that is more efficient, DW, 07/23/2018
# weight endpts with AK, DW, 07/23/2018
# do not weight endpts with AK, errors will be weighted by AK in cal.trans.err(), 
#   DW, 01/29/2019

endpts.trajfoot <- function(trajdat, timestr, ctmole.path){

  #### NOW grab a CO2.bio fluxes for each selected trajec ###
  # match 3-houly footprint and 3-houly biosperhic fluxes
  # first locate the date and hour for each backwards hours of STILT footprint
  recp.date <- as.POSIXlt(as.character(timestr), format = '%Y%m%d%H', tz = 'UTC')
  edp.time  <- tapply(trajdat$time, trajdat$indx, min)  # in mins

  # NOW, grab the endtraj for those min.time, using STRING MATCHING METHOD
  endtraj <- trajdat %>% group_by(indx) %>% filter(time == min(time)) %>%
                         mutate(date = recp.date + time * 60, 
                                timestr = substr(date, 1, 10),
                                yr  = as.numeric(substr(date, 1, 4)),
                                mon = as.numeric(substr(date, 6, 7)),
                                day = as.numeric(substr(date, 9, 10)),
                                hr  = as.numeric(substr(date, 12, 13))) %>% 
                         ungroup()
  uni.date <- unique(endtraj$timestr)

  # grab CT-NRT files for all backwards hours, and then put in a 3D fluxes array
  tmp <- nc_open(file.path(ctmole.path, list.files(ctmole.path, uni.date[1])))

  # grab variables
  ct.lat   <- ncvar_get(tmp, 'latitude') - 1
  ct.lon   <- ncvar_get(tmp, 'longitude') - 1.5
  ct.time  <- ncvar_get(tmp, 'time_components')	# UTC time components
  ct.level <- ncvar_get(tmp, 'level')	          # 25 levels
  ct.pres  <- ncvar_get(tmp, 'pressure')/100    # 26 boundaries for pressure
  ct.bound <- ncvar_get(tmp, 'boundary')	      # 26 boundaries
  rownames(ct.time) <- c('year', 'month', 'day', 'hour', 'minute', 'second')
  ct.hr <- ct.time['hour',] - 1 # every 3 hours

  # then match time to 3 hourly ct and get ct timestr for each particle,
  # find the correct hr and date.hr string given 3 hourly ct
  endtraj <- endtraj %>% mutate(match.hr = ct.hr[findInterval(hr, ct.hr)],
                                match.date.hr = paste(timestr, 
                                                      formatC(match.hr, 
                                                              width = 2, 
                                                              flag = 0)))

  # create a big 5D array for storing CO2 concentration, [lon, lat, 25layer, day, hour]
  ctco2.all <- array(0, dim = c(length(ct.lon), length(ct.lat), length(ct.level), 
                                length(uni.date), length(ct.hr)),
                        dimnames = list(ct.lon, ct.lat, ct.level, uni.date, ct.hr))
  
  ctgph.all <- array(0, dim = c(length(ct.lon), length(ct.lat), length(ct.bound), 
                                length(uni.date), length(ct.hr)),
                        dimnames = list(ct.lon, ct.lat, ct.bound, uni.date, ct.hr))

  for (f in 1 : length(uni.date)) {

      ctdat <- nc_open(list.files(ctmole.path, uni.date[f], full.names = T))

      # grab CO2 fields, 25 LEVELS FOR CO2
      ct.co2 <- ncvar_get(ctdat, 'co2')	# [LON, LAT, LEVEL, TIME]
      dimnames(ct.co2) <- list(ct.lon, ct.lat, ct.level, ct.hr)

      # grab geopotentail height, 26 BOUNDS FOR HGT
      ct.gph <- ncvar_get(ctdat, 'gph')	# [LON, LAT, BOUND, TIME]
      dimnames(ct.gph) <- list(ct.lon, ct.lat, ct.bound, ct.hr)

      # put into the big array for storing
      ctco2.all[,,,f,] <- ct.co2  # ppm
      ctgph.all[,,,f,] <- ct.gph  # in meter
      nc_close(ctdat)
  }  # end for f

  melt.co2 <- melt(ctco2.all)
  melt.gph <- melt(ctgph.all)
  colnames(melt.co2) <- list('lon', 'lat', 'layer', 'date', 'hr', 'co2')
  colnames(melt.gph) <- list('lon', 'lat', 'level', 'date', 'hr', 'gph')

  # to save time, subset CO2 concentration and gph, add buff of 2 deg
  melt.co2 <- melt.co2 %>% filter(lat >= floor(min(endtraj$lati) - 2), 
                                  lat <= ceiling(max(endtraj$lati) + 2),
                                  lon >= floor(min(endtraj$long) - 2), 
                                  lon <= ceiling(max(endtraj$long) + 2))

  melt.gph <- melt.gph %>% filter(lat >= floor(min(endtraj$lati) - 2), 
                                  lat <= ceiling(max(endtraj$lati) + 2), 
                                  lon >= floor(min(endtraj$long) - 2), 
                                  lon <= ceiling(max(endtraj$long) + 2))

  # after reading co2 mole fractions, assign lon.index, lat.index, day.index
  # and hour.index for each selected trajec
  sel.ct.lat <- unique(melt.co2$lat)
  sel.ct.lon <- unique(melt.co2$lon)

  endtraj <- endtraj %>% mutate(match.lati = sel.ct.lat[findInterval(lati, sel.ct.lat)],
                                match.long = sel.ct.lon[findInterval(long, sel.ct.lon)],
                                match.date = timestr, co2 = NA, 
                                zasl = zagl + zsfc)

  for (e in 1 : nrow(endtraj)) {

    # get CO2, based on matched lati, long, date, hr from 'endtraj'
    sel.co2 <- melt.co2 %>% filter(lat == endtraj$match.lati[e], 
                                   lon == endtraj$match.long[e], 
                                   as.character(date) == endtraj$match.date[e], 
                                   hr == endtraj$match.hr[e])

    sel.gph <- melt.gph %>% filter(lat == endtraj$match.lati[e], 
                                   lon == endtraj$match.long[e],
                                   as.character(date) == endtraj$match.date[e], 
                                   hr == endtraj$match.hr[e])

    # to further locate one CO2, requires ZASL from endtraj and melt.gph
    hgt.indx <- findInterval(endtraj$zasl[e], sel.gph$gph)

    # if below the lowest CT levels, use the first level
    if (hgt.indx == 0) hgt.indx <- 1

    # assign CO2 concentration from CT-NRT to endtraj
    endtraj$co2[e] <- (sel.co2 %>% filter(layer == hgt.indx))$co2
  }

  # weight endpts with AK, DW 07/23/2018
  uni.xhgt <- unique(endtraj$xhgt)
  endtraj <- endtraj %>% mutate(level = findInterval(xhgt, uni.xhgt))
  edp.co2 <- data.frame(indx = endtraj$indx, edp = endtraj$co2)

  return(edp.co2)
} # end of subroutine
