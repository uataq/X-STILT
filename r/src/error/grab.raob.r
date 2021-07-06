# readin NOAA radiosonde data and grab u, v wind speeds and wind directions
# by DW, 05/25/2017
# modify original script as a function, 08/29/2018

grab.raob <- function(raob.path, timestr, err.path, nhrs, 
                      format = c('ncdf', 'fsl')[1], overwrite = F){

  filename <- file.path(err.path, paste0(site, '_rad_', timestr, '.txt'))
  
  if ( file.exists(filename) & overwrite == F ) {
    options(scipen = 999)
    merge.var <- read.table(filename, sep = ',', header = T, stringsAsFactors = F)

  } else {

    # ------------------------- if netcdf format ------------------------------
    if (format == 'ncdf'){
      library(ncdf4)
      # grab the correct ncdf file
      raob.files <- list.files(path = raob.path, pattern = 'ncdf')
      raob.file  <- raob.files[grep(substr(timestr, 1, 8), raob.files)]

      # select times
      rel.time <- ncvar_get(dat, 'relTime')	# Seconds since (1970-1-1 00:00:0.0)
      rel.date <- as.POSIXct(rel.time, origin = '1970-01-01', tz = 'UTC')
      rel.timestr <- format(rel.date, '%Y%m%d%H%M%S')

      # grab station lat, lon, hgt, and name, and get rid of missing data
      site.name  <- ncvar_get(dat, 'staName')
      site.lat   <- ncvar_get(dat, 'staLat')
      site.lon   <- ncvar_get(dat, 'staLon')
      site.elev  <- ncvar_get(dat, 'staElev')
      merge.raob <- data.frame(site = site.name, lon = site.lon, lat = site.lat,
                               elev = site.elev, timestr = rel.timestr)
      merge.raob$nobs <- seq(1, nrow(merge.raob))

      # grabbing wind speeds, directions, pressure and sounding releasing times
      # all have dimensions [22 levels, total sites]
      grab.raob.var <- function(dat, var.id, varname){
        x <- ncvar_get(dat, var.id)
        melt.x <- melt(x)
        colnames(melt.x) <- list('nlev', 'nobs', varname)
        return(melt.x)
      }

      wind.mag <- grab.raob.var(dat, 'wsMan', 'ws')
      wind.dir <- grab.raob.var(dat, 'wdMan', 'wd')

      # mandatory pressure levels include sfc pres,
      # 1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20 mb
      pres  <- grab.raob.var(dat, 'prMan', 'pres')
      hgt   <- grab.raob.var(dat, 'htMan', 'hgt')
      temp  <- grab.raob.var(dat, 'tpMan', 'temp')

      # merge all variables
      merge.var <- Reduce(
        function(x, y) merge(x, y, by = c('nlev', 'nobs'), all = T),
        list(wind.mag, wind.dir, pres, hgt, temp)) %>%
        left_join(merge.raob, by = 'nobs') %>%
        filter(ws < 1E36, wd < 1E36, pres < 1E36, hgt < 1E36, temp < 1E36,
               lat <= 90, lon <= 180) %>%
        mutate(temp = temp - 273.15)
    } # end of ncdf format

    # ------------------------- if ASCII/FSL format ---------------------------
    if (format == 'fsl') {

      # grab the correct tmp file
      raob.files <- list.files(path = raob.path, pattern = 'tmp')
      raob.file  <- raob.files[grep(substr(timestr, 1, 8), raob.files)]

      # reading line by line, originate from JCL's
      dat  <- scan(file.path(raob.path, raob.file), sep = '\n', what = '')
      type <- as.numeric(substr(dat, 5, 7))   #type of identification line

      # grab time
      yr  <- as.numeric(substr(dat[type == 254], 35, 39))
      hr  <- as.numeric(substr(dat[type == 254], 13, 14))
      day <- as.numeric(substr(dat[type == 254], 20, 21))
      MON <- substr(dat[type == 254], 28, 30) # in character, change to number

      # use as.POSIXlt instead to simplify the code
      rel.date <- as.POSIXlt(paste(MON, day, yr, hr), format = '%b%d%Y%H')
      rel.timestr <- format(rel.date, '%Y%m%d%H%M%S')

      # grab locations and elev
      site.lat <- as.numeric(substr(dat[type == 1], 24, 28))
      xdir <- substr(dat[type == 1], 29, 29)
      site.lat[xdir == 'S'] <- -1 * site.lat[xdir == 'S']

      site.lon <- as.numeric(substr(dat[type == 1], 30, 35))
      xdir <- substr(dat[type == 1], 36, 36)
      site.lon[xdir == 'W'] <- -1 * site.lon[xdir == 'W']

      site.elev <- as.numeric(substr(dat[type == 1], 37, 42))
      merge.raob <- data.frame(lon = site.lon, lat = site.lat, elev = site.elev,
                               timestr = rel.timestr)

      # start to grab wind info
      # units of wind speed 'ms' = tenths of meters/sec; 'kt' = knots
      units <- substr(dat[type == 3], 48, 49)

      # start reading soundings
      SEL <- type == 254   #new sounding
      inds <- (1:length(dat))[SEL]
      inds <- c(inds, length(dat)+1)  # line index for a new sounding

      merge.var <- NULL
      for (i in 1:(length(inds) - 1)) {      #extract sounding, one at a time

        sel.dat <- dat[inds[i]:(inds[i+1] - 1)]

        #select only the data lines, 4 = mandatory level; 5 = significant level;
        # 6 = wind level (PPBB) (GTS or merged data);
        # 7 = tropopause level (GTS or merged data);
        # 8 = maximum wind level (GTS or merged data);
        # 9 = surface level
        sel <- !is.na(match(as.numeric(substr(sel.dat, 5, 7)),
                            c(9, 4, 5, 6, 7, 8)))
        dat.obs <- sel.dat[sel]

        # grab pressure, altitude, temp, dew points, wind dir and speeds
        pres <- as.numeric(substr(dat.obs, 10, 14)) * 0.1  # convert to mb
        hgt  <- as.numeric(substr(dat.obs, 17, 21))
        temp <- as.numeric(substr(dat.obs, 24, 28)) * 0.1
        wd   <- as.numeric(substr(dat.obs, 38, 42))
        ws   <- as.numeric(substr(dat.obs, 45, 49))

        # unit conversion
        if (units[i] == 'kt') ws <- ws * 1.61 * 1000 / (60 * 60)
        if (units[i] == 'ms') ws <- ws / 10
        #if unit is 'ms', then denotes TENTHS of [m/s]

        result <- data.frame(timestr = rel.timestr[i], lat = site.lat[i], 
                             lon = site.lon[i], elev = site.elev[i],
                             pres = pres, hgt = hgt, temp = temp, 
                             ws = ws, wd = wd)
        merge.var <- rbind(merge.var, result)
      }  # end for i

      # remove missing data
      merge.var <- merge.var %>% filter(temp < 1E3, ws < 1E3, wd < 1E4)
    }  # end if tmp

    # !!! degree true: from the true North to the wind barb,
    # so the actual wind direction is degree-true minus 180deg
    merge.var <- merge.var %>% mutate(u = sin((wd - 180) * pi/180) * ws,
                                      v = cos((wd - 180) * pi/180) * ws)

    # store in txt file
    write.table(merge.var, file = filename, sep = ',', quote = F, row.names = F)
  }  # end if file.exist

  # select nhrs back met
  date1 <- as.character(timestr)
  date2 <- strptime(date1, '%Y%m%d%H', tz = 'UTC')
  date2 <- format(date2 + nhrs * 60 * 60, '%Y%m%d%H')

  merge.var$timestr <- substr(merge.var$timestr, 1, 10)
  options(scipen = 0)
  sel.var <- merge.var %>% filter(timestr >= min(date1, date2), 
                                  timestr <= max(date1, date2))
  return(sel.var)
}

# end of script
