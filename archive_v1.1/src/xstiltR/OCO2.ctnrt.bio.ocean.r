# subroutine to readin CarbonTracker-NeatRealTime for OCO-2 project
# use CO2 optmized fluxes to get the biosperhic co2 exchange
# match weighted footprint with biosperhic co2 exchange
# input variables needed:
#       ident, footprint matrix, outpath for storing biospheric contributions
# Written by Dien Wu, 02/10/2017

# Update for CT-NRTv2017, DW
# combine bio with oceanic fluxes, DW, 05/03/2018

ctnrt.bio.ocean <- function(ident, foot, recp.info, ct.version, ctpath,
                            storeTF=FALSE, ncdfpath){

  library(ncdf4)

  # match 3-houly footprint and 3-houly biosperhic fluxes
  # first locate the date and hour for each backwards hours of STILT footprint
  foot.hour<- -1 * as.numeric(dimnames(foot)[[3]])-1

  # use weekdayhr() to return the acutal date,
  # weekdayhr<-function(yr,mon,day,hr,runtt,diffGMT=NA)
  # runtt needs to be in minutes
  back.times <- weekdayhr(yr=recp.info$recp.year, mon=recp.info$recp.mon,
                          day=recp.info$recp.day, hr=recp.info$recp.hour,
                          runtt=foot.hour*60)
  back.times <- data.frame(back.times)

  # fix actual hours, if back.hours%%!=0, 08/29/2017, DW
  if(recp.info$recp.hour %% 3 != 0){
    stand.hours <- seq(0, 21, 3)  # standard hour from CT, 00, 03,..., 21 UTC
    hr.index <- back.times$hr %% 3 != 0
    find.index <- findInterval(back.times[hr.index, "hr"], stand.hours)
    back.times[hr.index, "hr"] <- stand.hours[find.index]
  }

  # compute hours and days given footprint
  foot.hh <- paste(back.times$yr, formatC(back.times$mon, width=2, flag=0),
                                  formatC(back.times$day, width=2, flag=0),
                                  formatC(back.times$hr, width=2,flag=0),sep="")

  foot.dd <- paste(back.times$yr, formatC(back.times$mon, width=2, flag=0),
                              formatC(back.times$day, width=2, flag=0), sep="")
  foot.dd <- unique(foot.dd)

  # if foot hours differ from CT hours, the subroutine cannot find the CT times
  # break...
  if(unique(back.times$hr %% 3 != 0)){
    cat("ctnrt.bio.ocean(): foot hours differ from CT hours, NO files found...\n")

  }else{  # if not, we can find CT files

    # grab CT-NRT files at all backwards hours, and then put in 3D array
    # initialize 3D array for storing bio and oceanic fluxes,
    # only covers the regions in footprint [lat, lon, hours]

    # all lower left corner for footprint
    foot.lon <- as.numeric(dimnames(foot)[[2]])
    foot.lat <- as.numeric(dimnames(foot)[[1]])

    # use the same dimension as CT-NRT, NOT footprint,
    # remember to flip footprint dimension later
    ctbio <- array(0, dim=c(length(foot.lon), length(foot.lat),length(foot.hh)),
                      dimnames=list(foot.lon, foot.lat, foot.hh))
    ctocean <- ctbio
    store.length <- rep(0, length(unique(foot.dd)))

    # loop over each unique day and open CT files
    for (f in 1:length(foot.dd)){

      tmp.foot.hh <- foot.hh[substr(foot.hh,1,8) == foot.dd[f]]

      # open the daily file just once
      ctpattern <- paste("flux1x1.", foot.dd[f], sep="")
      ctfile <- list.files(path=ctpath, pattern=ctpattern)

      if(length(ctfile)==0){
        stop("ctnrt.bio.ocean(): NO CT flux files found...Please check...\n")
      }else{
        ctdat <- nc_open(file.path(ctpath, ctfile))
      }

      # 2D map for each 3 hours
      # 1 deg in N-S, 1 deg in E-W, move centered lat/lon to lower left
      # UTC time components
      if(ct.version == "v2016-1"){
        ct.lat <- ncvar_get(ctdat, "lat") - 0.5
        ct.lon <- ncvar_get(ctdat, "lon") - 0.5
        ct.time <- ncvar_get(ctdat, "date_components")
      } # if v2016

      if(ct.version == "v2017"){
        ct.lat <- ncvar_get(ctdat, "latitude") - 0.5
        ct.lon <- ncvar_get(ctdat, "longitude") - 0.5
        ct.time <- ncvar_get(ctdat, "time_components")
      } # if v2017

      rownames(ct.time) <- c("year", "mon", "day", "hour", "min", "sec")

      # convert to YYYYMMDDHHmm, move centered hours to the beginning hours
      # 1:30 --> 00:00 to 03:00, use the begining of each 3 hours
      ct.YYYYMMDDHHmmss <- paste(ct.time["year",],
                                 formatC(ct.time["mon",], width=2, flag=0),
                                 formatC(ct.time["day",], width=2, flag=0),
                                 formatC(ct.time["hour",]-1, width=2, flag=0),
                                 formatC(ct.time["min",]-30, width=2, flag=0),
                                 formatC(ct.time["sec",], width=2, flag=0),
                                 sep="")

      ct.hh <- substr(ct.YYYYMMDDHHmmss, 1, 10)

      # match footprint hr string with ct hr string, using match for hour index
      # also, select the lat, lon used in footprint, sel.ct.bio[lon, lat, hour]
      match.hh  <- match(tmp.foot.hh, ct.hh)
      match.lat <- match(foot.lat, ct.lat)
      match.lon <- match(foot.lon, ct.lon)

      # grab all biosperhic fluxes
      ct.bio <- ncvar_get(ctdat, "bio_flux_opt")  # unit in mol/m2/s, avg fluxes
      dimnames(ct.bio) <- list(ct.lon, ct.lat, ct.hh)
      ct.bio <- ct.bio * 1E6	# convert to umol/m2/s

      # grab all oceanic fluxes
      ct.ocean <- ncvar_get(ctdat, "ocn_flux_opt")	# unit in mol/m2/s
      dimnames(ct.ocean) <- list(ct.lon, ct.lat, ct.hh)
      ct.ocean <- ct.ocean * 1E6	# convert to umol/m2/s

      # select bio/ocean for certain lat, lon and time
      sel.ct.bio <- ct.bio[match.lon, match.lat, match.hh]
      sel.ct.ocean <- ct.ocean[match.lon, match.lat, match.hh]

      # create indices for storing
      store.length[f] <- length(tmp.foot.hh)  # store current files number
      if(f==1){
        min.store <- f
        max.store <- store.length[f]
      }else{
        min.store <- sum(store.length[1:(f-1)]) + 1
        max.store <- sum(store.length)
      }
      #print(c(min.store,max.store))

      # put into the big array for storing
      ctbio[,, seq(min.store, max.store)] <- sel.ct.bio  # umol/m2/s
      ctocean[,, seq(min.store, max.store)] <- sel.ct.ocean  # umol/m2/s
      nc_close(ctdat)
    }# end for f

    # Finally, we can match footprint with fluxes
    # flip footprint array first
    flip.foot <- aperm(foot, c(2,1,3))  # now [lon, lat, hr]
    dxco2.bio <- sum(flip.foot * ctbio)
    dxco2.ocean <- sum(flip.foot * ctocean)

    # store CT-NRT fluxes * 3D foot = CO2 contribution map into .nc file
    if(storeTF){

      ident2 <- gsub("&", "+", ident)
      cat("ctnrt.bio.ocean(): Storing foot x fluxes into ncdf files...\n")
      netcdf.name.bio <- paste("foot_bio_", ident2, ".nc", sep="")
      netcdf.name.ocean <- paste("foot_ocean_", ident2, ".nc", sep="")

      # foot.anthro x fluxes have dims of [LON, LAT]
      # then flip back to LAT, LON
      dco2.bio <- t(apply(flip.foot * ctbio, c(1,2), sum))
      dco2.ocean <- t(apply(flip.foot * ctocean, c(1,2), sum))

      contri.lat <- as.numeric(rownames(dco2.bio))
      contri.lon <- as.numeric(colnames(dco2.bio))

      #Set equal to our lat lon vectors we created earlier
      x <- ncdim_def("Lon", "degreesE", contri.lon)
      y <- ncdim_def("Lat", "degreesN", contri.lat)

      # flip 2D foot and store footprint in [LAT, LON]
      contri.var.bio <- ncvar_def(name="foot_bio", units="PPM", list(y,x),
                              longname="XCO2 change due to biospheric exchange")
      contri.var.ocean <- ncvar_def(name="foot_ocean", units="PPM", list(y,x),
                                 longname="XCO2 change due to oceanic exchange")

      # create nc files
      ncnew.bio <- nc_create(filename=netcdf.name.bio, vars=contri.var.bio)
      ncnew.ocean <- nc_create(filename=netcdf.name.ocean, vars=contri.var.ocean)

      #puts our variable into our netcdf file
      ncvar_put(nc=ncnew.bio, varid=contri.var.bio, vals=dco2.bio)
      ncvar_put(nc=ncnew.ocean, varid=contri.var.ocean, vals=dco2.ocean)

      nc_close(ncnew.bio)  #Closes our netcdf4 file
      nc_close(ncnew.ocean)

      #Move the output file name to our model output directory
      system(paste("mv", netcdf.name.bio, ncdfpath))
      system(paste("mv", netcdf.name.ocean, ncdfpath))
    }  # end store nc file

    return(c(dxco2.bio, dxco2.ocean))
  } # end if, checking whether ct files can be found

} # end of subroutine
