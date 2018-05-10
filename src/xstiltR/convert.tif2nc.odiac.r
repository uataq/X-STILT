# subroutine to readin tiff format of 1kmx1km ODIAC emissions
# and return nc format with selected emissions given a lat/lon domain
# DW, update with ODIACv2017, 11/01/2017

# add missing power plant in CARVA, PP10 and PP14 in ODIAC, DW, 11/28/2017

#zipTF  whether to unzip .tif.gz files
#addTF  whether to add the missing PP in CARVA

convert.tid2nc.odiac <- function(track.timestr, odiac.vname="2017", odiac.path,
                                 addTF=F, minlat, maxlat, minlon, maxlon){

  library(rgdal)
  library(sp)
  library(ncdf4)
  library(Hmisc)

  # get time
  YYYY <- substr(track.timestr, 1, 4)
  YYMM <- substr(track.timestr, 3, 6)
  YYYYMM <- substr(track.timestr, 1, 6)
  mod <- monthDays(as.Date(paste(YYYY, "-", substr(track.timestr, 5, 6),"-",
                                 substr(track.timestr, 7, 8), sep="")))

  # path and filename of the 1km GeoTiff ODIAC
  tiff.path <- file.path(odiac.path, "ODIAC_tiff", YYYY)
  gzfile <- list.files(path=tiff.path, pattern=YYMM)

  # check if RData file exist, no need to unzip and convert from tiff
  rd.file <- paste("ODIACv", odiac.vname, "_", YYYYMM, ".RData", sep="")
  rd <- file.path(odiac.path, rd.file)
  zipTF <- !file.exists(rd)

  if(zipTF){

    cat(paste("convert.tif2nc.odiac(): working on file", gzfile,"...\n"))
    if(grepl(".gz", gzfile)==TRUE){      # before reading tif file, unzip it

      # always unzip to odiac directory, do not unzip to group_data
      cat("Unzipping up...\n")
      system(paste("gunzip ", file.path(tiff.path, gzfile), sep=""))
      tiff.file <- substr(gzfile, 1, nchar(gzfile)-3)
    }else{
      tiff.file <- gzfile
    } # end if gz

    # Reading original tif data and convert to .RData
    #x <- GDAL.open(file.path(tiff.path, tiff.file))
    #info<-displayDataset(x)  # for plotting

    # after reading using readGDAL, the default dimension is (y,x)
    # 21600 rows and 43200 columns
    cat("Reading tiff file...\n")
    emiss.dat <- readGDAL(file.path(tiff.path, tiff.file))

    # assign as RData file
    cat("Assigning...\n")
    rd.file <- paste("ODIACv", odiac.vname, "_", YYYYMM, sep="")
    assignr(xname=rd.file, value=emiss.dat, path=odiac.path)

    # finally zip up the tif file
    cat("Zipping up...\n")
    system(paste("gzip ", file.path(tiff.path, tiff.file), "...\n"))
    gc()
  } # end if zipTF

  rd.file <- paste("ODIACv", odiac.vname, "_", YYYYMM, sep="")
  dat <- getr(path=odiac.path, xname=rd.file)

  # Centered Lat Lon, global, with 30 arcseond = 1/120 deg = 0.00833 deg res
  # 60 arcminutes in 1 degree; 60 arcseconds in 1 arcminute;
  # 1 arcsecond = 1/3600 degree !!!
  grid <- dat@grid
  lon.res <- grid@cellsize[1]
  lat.res <- grid@cellsize[2]
  minlon.ll <- grid@cellcentre.offset[1]-lon.res/2
  minlat.ll <- grid@cellcentre.offset[2]-lat.res/2
  maxlon.ll <- (grid@cells.dim[1]-1)*lon.res + minlon.ll
  maxlat.ll <- (grid@cells.dim[2]-1)*lon.res + minlat.ll

  # return the area map as [LON, LAT]
  cat("convert.tid2nc.odiac(): calculating the area...\n")
  emiss.lat <- seq(minlat, maxlat-lat.res, lat.res)
  emiss.lon <- seq(minlon, maxlon-lon.res, lon.res)
  area.co2 <- area(resolution=lon.res, start.lon=min(emiss.lon),
                   start.lat=min(emiss.lat), end.lon=max(emiss.lon),
                   end.lat=max(emiss.lat))
  t.area.co2 <- t(area.co2)

  # create a vector for lower left LAT LON
  lon.ll <- seq(minlon.ll, maxlon.ll, lon.res)
  lat.ll <- seq(minlat.ll, maxlat.ll, lat.res)

  # read CO2 data accordint to row first
  cat("reforming...\n")
  co2 <- matrix(dat@data$band1, nrow=21600, byrow=TRUE)

  # flip latitude, towards increasing trend,
  # since GDAL reads from UPPER LEFT CORNER, decreasing LAT
  fco2 <- co2[length(lat.ll):1, ]
  dimnames(fco2) <- list(lat.ll, lon.ll)

  # grab the Middle East region (Lat 0-50N, LON 0-60E) and put in .nc format
  # for PRD region (Lat 10-50N, LON 60-130E), 4800 y-pix, 8400 x-pix
  SEL.lat <- lat.ll >= minlat & lat.ll < maxlat
  SEL.lon <- lon.ll >= minlon & lon.ll < maxlon
  sel.co2 <- fco2[SEL.lat, SEL.lon]
  sel.lat <- as.numeric(dimnames(sel.co2)[[1]])
  sel.lon <- as.numeric(dimnames(sel.co2)[[2]])

  if(F){
    library(fields)
    log.sel.co2<-log10(sel.co2);log.sel.co2[log.sel.co2==-Inf]<-NA
    image.plot(sel.lon,sel.lat, t(log.sel.co2))
  }

  # convert the unit of CO2 emiss from Tonne Carbon/cell/month to umol/m2/s
  sel.co2 <- sel.co2*1E6/12*1E6 	# convert tonne-C to uomol-C (= umole-CO2)
  sel.co2 <- sel.co2/mod/24/60/60	# convert per month to per second
  sel.co2 <- sel.co2/t.area.co2		# convert per cell to per m2
  # NOW t.sel.co2 has unit of micromole-CO2/m2/s,
  # that can be used directly with footprint

  ### if addTF==TRUE, mannually add PP10 and PP14 in ODIAC
  if(addTF){
    pp.lat<-24.41796  # N
    pp.lon<-47.01769  # E, missing power plant to the SOuth of Riyadh

    # find the nearest ODIAC 1km*1km grid
    pp.lat.index <- findInterval(pp.lat, sel.lat)
    pp.lon.index <- findInterval(pp.lon, sel.lon)
    find.lat <- sel.lat[pp.lat.index] # pp will go to 24.41667N
    find.lon <- sel.lon[pp.lon.index] # pp will go to 47.01667E

    # original emission is only ~ 182 umol/m2/s
    find.co2 <- sel.co2[pp.lat.index, pp.lon.index]

    # get area for this grid box
    find.area <- t.area.co2[pp.lat.index, pp.lon.index]

    # riyadh.co2<-sel.co2[sel.lat>=22&sel.lat<=28,sel.lon>=45&sel.lon<=48]
    # no temporal sacling factor for now
    pp.co2 <- 3.5 # an assumed 3.5MtC/yr

    # convert Mt/yr to g/yr (1E12)
    # then to mol/yr and finally to umol-C/yr (=umol-CO2/yr)
    pp.co2 <- pp.co2*1E6*1E3*1E3/12*1E6
    pp.co2 <- pp.co2/365/24/60/60        # convert umol-CO2/yr to umol-CO2/s
    pp.co2 <- pp.co2/find.area           # convert umol-CO2/s to umol/m2/s

    sel.co2[pp.lat.index, pp.lon.index] <- find.co2 + pp.co2
  }

  #### storing as nc file
  # convert from .RData format to .nc format
  ll.string <- paste(minlat,"-",maxlat,"N_",minlon,"-",maxlon,"E", sep="")
  fn <- paste("odiac", odiac.vname, "_1kmx1km_",YYYYMM,"_", ll.string, sep="")

  netcdf.name <- paste(fn, ".nc",sep="")
  if(addTF)netcdf.name <- paste(fn, "_PP.nc",sep="")
  netcdf.name <- file.path(odiac.path, netcdf.name)

  # Set equal to our lat lon vectors we created earlier
  x <- ncdim_def( "lon",  "degreesE", sel.lon)
  y <- ncdim_def( "lat",  "degreesN", sel.lat)

  longname <- paste("Monthly mean CO2 ODIAC",odiac.vname," emission",sep="")
  if(addTF)longname <- paste("Monthly mean CO2 ODIAC",odiac.vname,
                            " emission with 2 power plants pp10+pp14", sep="")

  # store the emission grid in [LAT, LON]
  odiac.emiss<- ncvar_def(name="odiac_co2_emiss", units="micromol-co2/m2/s",
                          list(y,x),longname=longname)

  ncnew <- nc_create(filename=netcdf.name, vars=odiac.emiss)

  #puts our variable into our netcdf file
  ncvar_put(nc=ncnew, varid=odiac.emiss, vals=sel.co2)

  #Closes our netcdf4 file
  nc_close(ncnew)
  gc()

  return(netcdf.name)
}
