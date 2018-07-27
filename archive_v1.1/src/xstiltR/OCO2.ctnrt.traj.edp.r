#### subroutine to readin CarbonTracker-NeatRealTime for OCO-2 project
# use CO2_total mole fractions 3D fields to get the CO2 concentration
# for STILT particle endpoints --> trajec-endpoint boundary condition
# ALSO, NEED TO WEIGHT THE ENDPOINT CONCENTRATION USING AK & PWF !!!

# for current purpose (Riyadh, Middle East), only use global 3x2 files,
# may want to add a nested NAM 1x1 for US sites

# Written by Dien Wu, 09/12/2016

########################################################################
# updates --
# Before, we only consider CO2 boundary condition for model levels
# using trajec endpoints, so AK aloft is ZERO
# HOWEVER, what we should do is to use CT for upper atmos
#  + keep original OCO-2 AK aloft
# last fix on 04/06/2017, DW
#
# separate STILT level and CT-NRT levels, and output two XCO2, 04/24/2017, DW
# clean codes, DW, 06/12/2018
#########################################################################

# input variables needed--
#adjusted ak*pw profiles for STILT levels, .RData traj file, receptor time,

### list variable names for CO2_components,
# bg (background, component due to initial condition),
# bio (terrestrial biopshere exchange), ff (fossil fuel burning),
# fires (direct fire emissions), ocean (air-sea exchange),
# latitude, longitude, level, time (since 2000-01-01 00:00:00 UTC)

### list variable names for CO2_total,
# air_mass, blh (PBL thickness), co2, decimal_date,
# gph (geoponential_height), orography (sfc geopotential), pressure,
# specific_humidity, temperature, time_components

ctnrt.traj.edp <- function(ident, trajdat, recp.info, combine.profile,
                           ct.version, ctpath){

  # --------- PART 1: DEALING WITH ENDPOINTS OF TRAJ for model levels -------- #
  # for grabbing the CO2 concentration at endpoints,
  # we need the lon, lat, asl height and absolute time of STILT particles
  cat('ctnrt.traj.edp(): Grabbing trajec endpoints...It takes time...\n')

  # grab the information for particle that dropped off for each index,
  # (before some of them drop off when crossing ZERO lon)
  # change in 09/27/2016, following John's advice
  # find out the min time according to index number
  # as *ENDPOINTS* before DROPPING OFF, using TAPPLY commend
  min.time <- tapply(trajdat[, 'time'], trajdat[, 'index'], min)
  # return the min time [mins] back with attributes name as index

  # NOW, grab the endtraj for those min.time, using STRING MATCHING METHOD
  traj.time.index <- paste(trajdat[, 'time'], trajdat[, 'index'])
  end.time.index  <- paste(min.time, attributes(min.time)$dimnames[[1]])

  # returns row number for each index
  row.index <- match(end.time.index, traj.time.index)

  # finally grab trajec endpoints
  endtraj <- trajdat[row.index, ]

  # checking the time...
  #diff.time <- unique(endtraj[,'time'] - as.numeric(min.time))
  # should be zero, meaning grabbing endpoint correctly

  # THUS, endtraj data frame should always have numpar rows
  # compute the ASL heights for STILT particles
  asl <- endtraj[, 'agl'] + endtraj[, 'grdht']
  endtraj <- cbind(endtraj, asl)

  # use CO2_total for STILT background
  # locate which file to use, based on the receptor time and nhrs backwards
  # compute the end yr, mon, day, hour for nhrs back to open a CT file
  endp.time <- weekdayhr(
    yr = recp.info$recp.year, mon = recp.info$recp.mon,
		day = recp.info$recp.day, hr = recp.info$recp.hour,
		runtt = endtraj[, 'time'])    # runtt in mins

  # get unique day, year, mon for endpoints, and order them
  uni.endp.day  <- sort(unique(endp.time[,'day']))	 # e.g., day 26, 27, 28
  uni.endp.year <- sort(unique(endp.time[,'yr']))
  uni.endp.mon  <- sort(unique(endp.time[,'mon']))

  # since CT files contain 3D CO2 fields every 3 hours in UTC,
  # assign a endpoint time index to each particle
  # CT time: centered time on 1:30, 4:30, 7:30, 10:30, 13:30, 16:30, 19:30 & 22:30
  # time range from 00-03, 03-06, 06-09, 09-12, 12-15, 15-18, 18-21 and 21-24.

  # Open all CT files into a list before assigning background CO2
  all.ct.gph <- rep(list(NULL), length(uni.endp.day))
  all.ct.co2 <- all.ct.gph

  for (f in 1:length(uni.endp.day)) {

    # find and open the CT file
    ctpattern <- paste0(
			formatC(uni.endp.year, width = 2, flag = 0), '-',
		  formatC(uni.endp.mon, width = 2, flag = 0), '-',
		  formatC(uni.endp.day[f], width = 2, flag = 0), '.nc')
    ctfile <- list.files(path = ctpath, pattern = ctpattern)

    cat(paste('ctnrt.traj.edp(): Reading from file', ctfile, '...\n'))
    ctdat <- nc_open(file.path(ctpath, ctfile))

    # ---------------- readin CT-NRT for getting molefractions --------------- #
    # 2 deg in N-S, 3 deg in E-W, move centered lat, lon to lower left
    ct.lat  <- ncvar_get(ctdat, 'latitude') -1
    ct.lon  <- ncvar_get(ctdat, 'longitude')-1.5
    ct.time <- ncvar_get(ctdat, 'time_components')	# UTC time components
    rownames(ct.time)<-c('year', 'mon', 'day', 'hour', 'min', 'sec')

    # convert to YYYYMMDDHHmmss
    ct.timestr <- paste0(
			ct.time['year',], formatC(ct.time['mon',],  width = 2, flag = 0),
			formatC(ct.time['day',],  width = 2, flag = 0),
			formatC(ct.time['hour',], width = 2, flag = 0),
			formatC(ct.time['min',],  width = 2, flag = 0),
			formatC(ct.time['sec',],  width = 2, flag = 0))

    ct.level <- ncvar_get(ctdat, 'level')	       # 25 levels
    ct.bound <- ncvar_get(ctdat, 'boundary')	   # 26 boundaries

    # grab geopotentail height, 26 BOUNDS FOR HGT
    ct.gph <- ncvar_get(ctdat, 'gph')	           # [LON, LAT, BOUND, TIME]
    dimnames(ct.gph) <- list(ct.lon, ct.lat, ct.bound, ct.timestr)

    # grab CO2 fields, 25 LEVELS FOR CO2
    ct.co2 <- ncvar_get(ctdat, 'co2')	           # [LON, LAT, LEVEL, TIME]
    dimnames(ct.co2) <- list(ct.lon, ct.lat, ct.level, ct.timestr)

    # combine all arrays into list
    all.ct.gph[[f]] <- ct.gph
    all.ct.co2[[f]] <- ct.co2

    nc_close(ctdat)

  } # end f for reading in ct files

  names(all.ct.gph) <- uni.endp.day
  names(all.ct.co2) <- uni.endp.day

  ## use lat/lon/asl to locate the CT grid (all CT files have same lat lon grid)
  lat.index <- findInterval(endtraj[, 'lat'], ct.lat)
  lon.index <- findInterval(endtraj[, 'lon'], ct.lon)

  # CHECKING... diff.lat < 2, diff.lon < 3
  #diff.lat <- endtraj[, 'lat'] - ct.lat[lat.index]
  #diff.lon <- endtraj[, 'lon'] - ct.lon[lon.index]

  # locate which day of this endpoint locates
  last.day <- endp.time[, 'day']

  # find when CT hour range that the last STILT hour (endp.time[,'hr']) falls
  # CT hours: 0130 (for 0000-0300), 0430, 0730, 1030, 1330, 1630, 1930, and 2230
  # calculate the CT hour index, and grab the CO2 fields at this hour
  last.hour <- endp.time[, 'hr']
  ct.hour.index <- trunc(last.hour / 3 + 1)	# which ct time to look for CO2

  # ---------- DEALING WITH CT 26 BOUNDARY AND 25 LEVEL for hgt.index -------- #
  # assign 'hgt.index' to each particles, referring which CT levels to grab,
  # by comparing ASL from STILT and GPH from CT
  # our goal is to relate CO2 concentration with the hgts, not levels
  # Determines the location, i.e., index of the (first) min/max of
  # a numeric (or logical) vector, using 'which.min' function
  cat('ctnrt.traj.edp(): Grabbing CT CO2 for each particle...\n')

  # initial ct.bg.co2 with NA
  ct.bg.co2 <- rep(NA, nrow(endtraj))
  endtraj   <- cbind(endtraj, ct.bg.co2)

  # loop over each endpoint
  for (p in 1:nrow(endtraj)) {

    # simplify on 04/12/2017, DW
    p.asl <- endtraj[p, 'asl']  # ASL height for particle p

    # locate the day first, all.ct.*-->ct.*.day
    ct.gph.day <- all.ct.gph[names(all.ct.gph) == last.day[p]][[1]]
    ct.co2.day <- all.ct.co2[names(all.ct.co2) == last.day[p]][[1]]

    # then locate an air column, based on lat/lon/hour, ct.*.day-->ct.*.loc
    ct.gph.loc <- ct.gph.day[lon.index[p], lat.index[p], , ct.hour.index[p]]
    ct.co2.loc <- ct.co2.day[lon.index[p], lat.index[p], , ct.hour.index[p]]

    # compute boundary index
    bound.index <- findInterval(p.asl, ct.gph.loc)

    # if bound.index==0, meaning particel hgt is below the first CT level,
    # simply assign bound.index to 1, use 1st CT level CO2
    if(bound.index == 0)bound.index <- 1

    # convert to level index
    # since level 1 CO2 is the concentration between BOUND 1 and BOUND 2
    # so, level index == lower boundary index
    # grab co2 from 'ct.co2.loc' according to 'bound.index'
    endtraj[p, 'ct.bg.co2'] <- ct.co2.loc[bound.index]
  } # end loop p assigning CO2 values

  # check: whether ct.bg.co2 all have values, if zero, good
  #check <- length(which(is.na(endtraj[,'ct.bg.co2']) == TRUE))

  # take the mean of ensemble at each STILT releasing levels
  co2.bg.stilt <- tapply(endtraj[, 'ct.bg.co2'], endtraj[, 'level'], mean)

  # because ct.bg.co2 is not normal distributed,
  # try 50th percentile, instead of mean, DW, 11/16/2017
  #co2.bg.stilt <- tapply(endtraj[, 'ct.bg.co2'], endtraj[, 'level'], median)

  # --------- PART 2: DEALING WITH ENDPOINTS OF TRAJ for model levels -------- #
  # grab CT pressure and interpolate CT at receptors for high up
  cat('ctnrt.traj.edp(): Grabbing CT CO2 above max model level...\n')

  # grab CT file right at the receptor
  ctpattern <- paste0(
		formatC(recp.info$recp.year, width = 2, flag = 0), '-',
    formatC(recp.info$recp.mon,  width = 2, flag = 0), '-',
    formatC(recp.info$recp.day,  width = 2, flag = 0), '.nc')

  ctfile <- list.files(path = ctpath, pattern = ctpattern)
  ctdat <- nc_open(file.path(ctpath, ctfile))

  # find the closest lat lon hour grid (to receptor) and grab the CT profiles
  # and then interpolate onto combine levels
  # 2 deg in N-S, 3 deg in E-W, move centered lat, lon to lower left
  ct.lat  <- ncvar_get(ctdat, 'latitude') -1
  ct.lon  <- ncvar_get(ctdat, 'longitude')-1.5
  ct.time <- ncvar_get(ctdat, 'time_components')	# UTC time components
  rownames(ct.time) <- c('year', 'mon', 'day', 'hour', 'min', 'sec')

  # move centered time to the beginning of each time period
  ct.time['hour', ] <- ct.time['hour', ] - 1   # minus 1 hr
  ct.time['min', ]  <- ct.time['min', ] - 30   # minus 30 mins
  ct.hour <- ct.time['hour', ]	               # starting hours

  # compare recptor hour with CT hour, and find the time index
  # similar for lat, lon
  hr.index  <- findInterval(recp.info$recp.hour, ct.hour)
  lat.index <- findInterval(recp.info$recp.lat, ct.lat)
  lon.index <- findInterval(recp.info$recp.lon, ct.lon)

  # grab CO2 fields, 25 LEVELS FOR CO2
  ct.co2   <- ncvar_get(ctdat, 'co2')	         # [LON, LAT, LEVEL, TIME]
  ct.level <- ncvar_get(ctdat, 'level')	       # 25 levels
  dimnames(ct.co2) <- list(ct.lon, ct.lat, ct.level, ct.hour)

  # grab pressure level instead of geopotentail height, at 26 boundaries
  # convert Pa, to hPa, [LON, LAT, BOUND, TIME]
  ct.pres <- ncvar_get(ctdat, 'pressure') / 100
  ct.bound <- ncvar_get(ctdat, 'boundary')	  # 26 boundaries
  dimnames(ct.pres) <- list(ct.lon, ct.lat, ct.bound, ct.hour)

  # NOW, grab the CT CO2 profiles at the recptor
  sel.co2  <- ct.co2[lon.index, lat.index, , hr.index]  # 25 CT levels
  sel.pres <- ct.pres[lon.index, lat.index, , hr.index]  # pressure in hPa

  # interpolating CT CO2 profiles onto combine levels, only above STILT levels
  model.pres <- combine.profile$pres
  max.model.pres <- model.pres[max(which(combine.profile$stiltTF == T))]
  upper.pres <- model.pres[model.pres < max.model.pres]

  # assigning CT CO2 values to each model pressure, go through each
  co2.bg.upper <- NULL
  for (l in 1:length(upper.pres)) {

	  # locate which layer that model level falls in
	  diff.pres <- sel.pres - upper.pres[l]

	  # differ from GPH, pressure is opposite
	  pres.bound <- as.numeric(attributes(diff.pres[diff.pres >= 0])$names)

    # --> lower boundary index, if pres.bound==0, use sfc level
    if (length(pres.bound) == 0) {
      bound.index <- 1
    } else{
      bound.index <- max(pres.bound)
    }  # end if
	  co2.bg.upper <- c(co2.bg.upper, sel.co2[names(sel.co2) == bound.index])

  } # end for l, assigning CO2 background for upper levels


  # ------------ PART 3: COMBINING UPPER AND LOWER bg CO2 ------------------- #
  co2.bg.combine <- c(co2.bg.stilt, co2.bg.upper)
  combine.profile$ctnrt.bg <- co2.bg.combine

  # apply ak pw weighting for each level mean
  xco2.bg <- sum(combine.profile$ak.pw * combine.profile$ctnrt.bg, na.rm=T)

  return(list(combine.profile,xco2.bg))
}

# end of subroutine