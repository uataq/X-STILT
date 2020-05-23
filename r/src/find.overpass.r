# subroutine to search OCO-2 files for any overpases for a given region
# DW, 05/15/2017
# add one count for good quality data, DW, 12/20/2017
# 'date.range' for date range, c(YYYYMMDD_1, YYYYMMDD_2)
# 'oco2.path' default in lin-group

# lon.lat contains c(minlat, maxlat, minlon, maxlon, ...), NEED THIS ORDER!!!
# default city center is for Riyadh
# change the order of 'lon.lat',
# --- c(min.lon, max.lon, min.lat, max.lat, city.lon, city.lat), DW
# add a flag ''urbanTF' for searching soundings near urban region, DW, 06/15/2018
# dlat, dlon for lat,lon from city center
# update for v9 data, DW, 10/19/2018 
# drop the scientific notation for sounding ID, DR, DW, 09/04/2019

find.overpass <- function(date.range, lon.lat, oco2.ver = c('b7rb','b8r', 'b9r')[3],
                          oco2.path, urbanTF = F, dlon = 0.5, dlat = 0.5){

  library(geosphere); library(ncdf4); library(dplyr)

  # path and filename for storing OCO-2 info
  all.file <- list.files(pattern = 'oco2_LtCO2_', path = oco2.path)

  # get rid of some characters
  file.info <- gsub('oco2_LtCO2_', '', all.file)
  file.info <- gsub('.nc4', '', file.info)
  file.info <- strsplit.to.df(file.info)
  #file.info <- data.frame(matrix(unlist(strsplit(file.info, '_')), ncol = 3,
  #                               byrow = T), stringsAsFactors = F)
  colnames(file.info)[1] <- 'timestr'
  all.timestr <- file.info$timestr

  # for b8 or b9
  if (oco2.ver == 'b7rb' | oco2.ver == 'b10') {
    all.timestr[nchar(all.timestr) == 6] <-
      paste0('20', all.timestr[nchar(all.timestr) == 6])

    oco2.file <- all.file[all.timestr >= date.range[1] & all.timestr <= date.range[2]]
    timestr <- all.timestr[all.timestr >= date.range[1] & all.timestr <= date.range[2]]

  } else if (oco2.ver != 'b7rb') {
    SEL.day <- all.timestr >= substr(date.range[1], 3, 8) &
               all.timestr <= substr(date.range[2], 3, 8)
    oco2.file <- all.file[SEL.day]
    timestr <- paste0('20', substr(oco2.file, 12, 17))
  } 

  # loop over each overpass
  result <- NULL
  for (f in 1:length(oco2.file)) {

    if (f %% 25 == 0)
    cat(paste('#-- ', signif(f / length(oco2.file), 3) * 100, '% SEARCHED --#\n'))
    dat <- nc_open(file.path(oco2.path, oco2.file[f]))

    ## grabbing OCO-2 levels, lat, lon
    oco2.lev <- ncvar_get(dat, 'levels')
    oco2.lat <- ncvar_get(dat, 'latitude')
    oco2.lon <- ncvar_get(dat, 'longitude')
    xco2 <- ncvar_get(dat, 'xco2'); xco2[xco2 == -999999] <- NA
    qf <- ncvar_get(dat, 'xco2_quality_flag')

    # drop the scientific notation for sounding ID, DR, DW, 09/04/2019
    id <- format(ncvar_get(dat, 'sounding_id'), scientific = F)
    obs <- data.frame(lat = as.numeric(oco2.lat), lon = as.numeric(oco2.lon), 
                      qf = as.numeric(qf), xco2 = as.numeric(xco2), 
                      timestr = as.numeric(substr(id, 1, 10)), 
                      stringsAsFactors = F) 
    
    # Warn level being removed for lite v9 data, DW, 10/15/2018 
    if (oco2.ver %in% c('b8r', 'b7rb')) {
      wl <- ncvar_get(dat, 'warn_level')
      obs <- obs %>% mutate(wl = wl)
    } # end if not v9

    obs <- obs %>% filter(lat >= lon.lat$minlat & lat <= lon.lat$maxlat &
                          lon >= lon.lat$minlon & lon <= lon.lat$maxlon) 

    # count overpasses
    tot.count <- nrow(obs)
    qf.count  <- nrow(obs %>% filter(qf == 0))
    if (oco2.ver == 'b8r') wl.count <- nrow(obs %>% filter(wl == 0))
    if (oco2.ver == 'b7rb') wl.count <- nrow(obs %>% filter(wl <= 15))

    nc_close(dat)
    
    # store results if there are soundings over
    if (tot.count > 0) {
      uni.timestr <- unique(obs$timestr)

      # also search for soundings near city center,
      if (urbanTF) {
        urban.dat <- obs %>% filter(lon >= (lon.lat$citylon - dlon),
                                    lon <= (lon.lat$citylon + dlon),
                                    lat >= (lon.lat$citylat - dlat),
                                    lat <= (lon.lat$citylat + dlat))
        tot.urban.count <- nrow(urban.dat)
        qf.urban.count  <- nrow(urban.dat %>% filter(qf == 0))
        if (oco2.ver == 'b8r') wl.urban.count <- nrow(urban.dat %>% filter(wl == 0))
        if (oco2.ver == 'b7rb') wl.urban.count <- nrow(urban.dat %>% filter(wl <= 15))

        # combine 
        tmp <- data.frame(timestr = as.numeric(uni.timestr), 
                          tot.count = as.numeric(tot.count),
                          qf.count = as.numeric(qf.count), 
                          tot.urban.count = as.numeric(tot.urban.count), 
                          qf.urban.count = as.numeric(qf.urban.count),
                          stringsAsFactors = F)
        if (oco2.ver %in% c('b8r', 'b7rb')) 
          tmp <- cbind(tmp, wl.count = as.numeric(wl.count), 
                            wl.urban.count = as.numeric(wl.urban.count))
      
      } else {
        tmp <- data.frame(timestr = as.numeric(uni.timestr), 
                          tot.count = as.numeric(tot.count),
                          qf.count = as.numeric(qf.count), stringsAsFactors = F)
        
        if (oco2.ver %in% c('b8r', 'b7rb')) 
          tmp <- cbind(tmp, wl.count = as.numeric(wl.count))
      } # end if urbanTF

      result <- rbind(result, tmp)
    } else next  # if no data, jump to next file
     # end if tot.count

  }  # end for f

  result <- result[order(result$timestr),]
  return(result)
}
