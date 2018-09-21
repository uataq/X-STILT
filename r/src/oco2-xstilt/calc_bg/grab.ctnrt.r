# subroutine to grab biospheric or oceanic fluxes from CT-NRT, 
# according to footprint backwards hours
# DW, 09/15/2018

# foot.dd for unique time string, YYYYMMDD, can be a vector
# foot.hh for all backward time strings, YYYYMMDDHH, a vector
# both above variables should have the same order of foot time dim

grab.ctnrt.flux <- function(foot, foot.lat, foot.lon, foot.dd, foot.hh, 
  flux.path, ct.ver) {

    # use the same dimension as CT-NRT, NOT footprint,
    # remember to flip footprint dimension later
    ctbio <- array(0, dimnames = dimnames(foot), dim = dim(foot))
    ctocn <- ctbio
    store.length <- rep(0, length(foot.dd))

    # loop over each unique day and open CT files
    for (f in 1 : length(foot.dd)){

      tmp.foot.hh <- foot.hh[substr(foot.hh, 1, 8) == foot.dd[f]]

      # open the daily file just once
      flux.file <- list.files(path = flux.path, 
        pattern = paste0('flux1x1.', foot.dd[f]), full.names = T)
      if (length(flux.file) == 0) {
        stop('ctnrt.bio.ocean(): NO CT flux files found...Please check...\n')
      } else {
        ctdat <- nc_open(flux.file)
      } # end if 

      # 2D map for every 3 hours, 1 deg in N-S, 1 deg in E-W
      # no need to convert to centered lat lon, as foot lat lon is centered
      # UTC time components      
      if (ct.ver == 'v2016-1') {
        ct.lat  <- ncvar_get(ctdat, 'lat')
        ct.lon  <- ncvar_get(ctdat, 'lon')
        ct.time <- ncvar_get(ctdat, 'date_components')
      } # if v2016

      if (ct.ver == 'v2017') {
        ct.lat  <- ncvar_get(ctdat, 'latitude')
        ct.lon  <- ncvar_get(ctdat, 'longitude')
        ct.time <- ncvar_get(ctdat, 'time_components')
      } # if v2017
      rownames(ct.time) <- c('year', 'mon', 'day', 'hour', 'min', 'sec')
      ct.time <- data.frame(t(ct.time))

      # convert to YYYYMMDDHHmm, move centered hours to the beginning hours
      # 1:30 --> 00:00 to 03:00, use the begining of each 3 hours
      ct.YYYYMMDDHHmmss <- paste0(ct.time$year,
                                 formatC(ct.time$mon, width = 2, flag = 0),
                                 formatC(ct.time$day, width = 2, flag = 0),
                                 formatC(ct.time$hour - 1, width = 2, flag = 0),
                                 formatC(ct.time$min - 30, width = 2, flag = 0),
                                 formatC(ct.time$sec, width = 2, flag = 0))
      ct.hh <- substr(ct.YYYYMMDDHHmmss, 1, 10)

      # match footprint hr string with ct hr string, using match for hour index
      # also, select the lat, lon used in footprint, sel.ct.bio[lon, lat, hour]
      match.hh  <- match(tmp.foot.hh, ct.hh)
      match.lat <- match(foot.lat, ct.lat) # centered lat, lon 
      match.lon <- match(foot.lon, ct.lon)

      # grab all biosperhic fluxes
      ct.bio <- ncvar_get(ctdat, 'bio_flux_opt')  # unit in mol/m2/s, avg fluxes
      dimnames(ct.bio) <- list(ct.lon, ct.lat, ct.hh)
      ct.bio <- ct.bio * 1E6	# convert to umol/m2/s

      # grab all oceanic fluxes
      ct.ocean <- ncvar_get(ctdat, 'ocn_flux_opt')	# unit in mol/m2/s
      dimnames(ct.ocean) <- list(ct.lon, ct.lat, ct.hh)
      ct.ocean <- ct.ocean * 1E6	# convert to umol/m2/s

      # select bio/ocean for certain lat, lon and time
      sel.ct.bio   <- ct.bio[match.lon, match.lat, match.hh]
      sel.ct.ocean <- ct.ocean[match.lon, match.lat, match.hh]

      # create indices for storing
      store.length[f] <- length(tmp.foot.hh)  # store current files number
      if (f == 1) {
        min.store <- f
        max.store <- store.length[f]
      } else {
        min.store <- sum(store.length[1 : (f - 1)]) + 1
        max.store <- sum(store.length)
      }
      #print(c(min.store,max.store))

      # put into the big array for storing
      ctbio[,, seq(min.store, max.store)] <- sel.ct.bio  # umol/m2/s
      ctocn[,, seq(min.store, max.store)] <- sel.ct.ocean  # umol/m2/s
      nc_close(ctdat)
    } # end for f

    ct.flux <- list(bio = ctbio, ocn = ctocn)
    return(ct.flux)
}

