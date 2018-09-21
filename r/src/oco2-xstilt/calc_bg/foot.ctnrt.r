#### subroutine to readin CT-NRT fluxes and then couple with STILT footprint
# as footprints have already been weighted by AK and PW,
# just multiple fluxion with 2D footprint map
# written by Dien Wu, 09/13/2016

foot.ctnrt <- function(foot.file, ct.ver, flux.path, timestr, nhrs, 
  txtfile = NULL, writeTF = F){

  library(ncdf4)

  # calculate backtime
  recp.time <- as.POSIXct(timestr, format = '%Y%m%d%H', tz = 'UTC')
  foot.time <- seq(sign(nhrs), nhrs, sign(nhrs))
  back.time <- recp.time + foot.time * 60 * 60  # nhrs with sign
  back.timestr <- format(back.time, format = '%Y%m%d%H')
  back.times <- data.frame(yr  = as.numeric(substr(back.timestr, 1, 4)), 
                           mon = as.numeric(substr(back.timestr, 5, 6)),
                           day = as.numeric(substr(back.timestr, 7, 8)),
                           hr  = as.numeric(substr(back.timestr, 9, 10)))

  # fix actual hours to CT hours, if back.hours %% !=0, 08/29/2017, DW
  if (as.numeric(substr(timestr, 9, 10)) %% 3 != 0){
    std.hr <- seq(0, 21, 3)  # standard hour from CT, 00, 03,..., 21 UTC
    hr.index <- back.times$hr %% 3 != 0
    find.index <- findInterval(back.times[hr.index, 'hr'], std.hr)
    back.times[hr.index, 'hr'] <- std.hr[find.index]
  }

  # compute hours and days given footprint
  foot.hh <- paste0(back.times$yr, 
                    formatC(back.times$mon, width = 2, flag = 0),
                    formatC(back.times$day, width = 2, flag = 0),
                    formatC(back.times$hr, width = 2, flag = 0))
  foot.dd <- unique(substr(foot.hh, 1, 8))

  # from foot.file, get receptor info
  receptor <- unlist(strsplit(gsub('_X_1x1_foot.nc', '', basename(foot.file)), '_'))
  receptor <- as.data.frame(matrix(receptor, byrow = T, ncol = 3),
    stringsAsFactors = F) %>% mutate_all(funs(as.numeric), colnames(receptor))
    # mutate_all() convert character to numberic
  colnames(receptor) <- list('timestr', 'lon', 'lat')

  order.index <- order(receptor$lat)
  receptor    <- receptor[order.index, ]
  foot.file   <- foot.file[order.index]
  receptor$xco2.bio <- NA
  receptor$xco2.ocn <- NA

  # read in footprint
  foot.dat <- nc_open(foot.file[1])
  foot <- ncvar_get(foot.dat, 'foot')
  foot.lon <- ncvar_get(foot.dat, 'lon')
  foot.lat <- ncvar_get(foot.dat, 'lat')
  dimnames(foot) <- list(foot.lon, foot.lat, ncvar_get(foot.dat, 'time')) 
  nc_close(foot.dat)

  # grab CT fluxes 
  # call function to grab biospheric or oceanic fluxes 
  #source('r/dependencies.r') # source all functions
  ct.flux <- grab.ctnrt.flux(foot, foot.lat, foot.lon, foot.dd, foot.hh, 
    flux.path, ct.ver)
  ct.bio  <- ct.flux$bio  # with centered lat, lon for fluxes
  ct.ocn  <- ct.flux$ocn 

  # then loop over each receptor
  for (r in 1:nrow(receptor)) {

    # read in footprint
    foot.dat <- nc_open(foot.file[r])
    foot     <- ncvar_get(foot.dat, 'foot')
    foot.lon <- ncvar_get(foot.dat, 'lon')
    foot.lat <- ncvar_get(foot.dat, 'lat')
    dimnames(foot) <- list(foot.lon, foot.lat, ncvar_get(foot.dat, 'time')) 

    # Finally, we can match footprint with fluxes
    if (class(all.equal(dim(foot), dim(ct.bio))) == 'logical') {
      xco2.bio <- foot * ct.bio
      xco2.ocn <- foot * ct.ocn
      receptor$xco2.bio[r] <- sum(xco2.bio)
      receptor$xco2.ocn[r] <- sum(xco2.ocn)
      #print(receptor$xco2.bio[r])

       # store into the same workding dir
      fn.bio <- file.path(dirname(foot.file[r]),
        gsub('1x1_foot', 'foot_bio', basename(foot.file[r])))
      fn.ocn <- file.path(dirname(foot.file[r]),
        gsub('1x1_foot', 'foot_ocn', basename(foot.file[r])))

      time_out <- ncvar_get(foot.dat, 'time')
      lon <- foot.lon
      lat <- foot.lat
      write_contri(xco2.bio, time_out, lon, lat, fn.bio, str = 'bio')
      write_contri(xco2.ocn, time_out, lon, lat, fn.ocn, str = 'ocn')
    } else {next} # else if dimension not match
   
    nc_close(foot.dat)
  }  # end for r

  # finally, write in a txt file
  if (writeTF)
    write.table(x = receptor, file = txtfile, sep = ',', row.names = F, quote = F)

  return(receptor)
} # end of subroutine
