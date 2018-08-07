# subroutine to search OCO-2 files for any overpases for a given region
# DW, 05/15/2017
# add one count for good quality data, DW, 12/20/2017
# 'date.range' for date range, c(YYYYMMDD_1, YYYYMMDD_2)
# 'oco2.path' default in lin-group

# target.region contains c(min.lat, max.lat, min.lon, max.lon), NEED THIS ORDER!!!
# default city center is for Riyadh
# change the order of 'target.region',
# --- c(min.lon, max.lon, min.lat, max.lat, city.lon, city.lat), DW
# add a flag ''urbanTF' for searching soundings near urban region, DW, 06/15/2018
# dlat, dlon for lat,lon from city center

find.overpass <- function(date.range, target.region,
  oco2.ver = c('b7rb','b8r'),
  oco2.path = file.path('/uufs/chpc.utah.edu/common/home/lin-group5/wde/input_data/OCO-2/L2',
              paste0('OCO2_lite_', oco2.ver)),
  urbanTF = F, dlon = 0.5, dlat = 0.5){

  library(geosphere); library(ncdf4); library(dplyr)

  # path and filename for storing OCO-2 info
  all.file <- list.files(pattern = 'oco2_LtCO2_', path = oco2.path)

  # get rid of some characters
  file.info <- gsub('oco2_LtCO2_', '', all.file)
  file.info <- gsub('.nc4', '', file.info)
  file.info <- data.frame(matrix(unlist(strsplit(file.info, '_')), ncol = 3,
    byrow = T), stringsAsFactors = F)
  colnames(file.info) <- c('timestr', 'ver', 'type')
  all.timestr <- file.info$timestr


  if (oco2.ver == 'b8r') {
    SEL.day <- all.timestr >= substr(date.range[1], 3, 8) &
               all.timestr <= substr(date.range[2], 3, 8)
    oco2.file <- all.file[SEL.day]
    timestr <- paste0('20', substr(oco2.file, 12, 17))

  } else if (oco2.ver == 'b7rb') {

    all.timestr[nchar(all.timestr) == 6] <-
      paste0('20', all.timestr[nchar(all.timestr) == 6])

    oco2.file <- all.file[all.timestr >= date.range[1] & all.timestr <= date.range[2]]
    timestr <- all.timestr[all.timestr >= date.range[1] & all.timestr <= date.range[2]]
  } # end if oco2.ver


  # loop over each overpass
  result <- NULL

  for (f in 1:length(oco2.file)) {

    if (f %% 25 == 0)
    cat(paste('#-- ', signif(f/length(oco2.file), 3)*100, '% SEARCHED --#\n'))

    oco2.dat <- nc_open(file.path(oco2.path, oco2.file[f]))
    oco2.lat <- ncvar_get(oco2.dat, 'latitude')
    oco2.lon <- ncvar_get(oco2.dat, 'longitude') # grabbing OCO-2 levels, lat, lon
    oco2.qf  <- ncvar_get(oco2.dat, 'xco2_quality_flag')
    oco2.wl  <- ncvar_get(oco2.dat, 'warn_level')

    # determine whether there are any overpasses for the target region
    SEL <- oco2.lon >= target.region[1] & oco2.lon <= target.region[2] &
           oco2.lat >= target.region[3] & oco2.lat <= target.region[4]
    tot.count <- length(which(SEL))
    qf.count  <- length(which(SEL & oco2.qf == 0))
    if (oco2.ver == 'b8r')  wl.count <- length(which(SEL & oco2.wl == 0))
    if (oco2.ver == 'b7rb') wl.count <- length(which(SEL & oco2.wl < 15))

    # also search for soundings near city center,
    # city.lon = target.region[5]; city.lat = target.region[6]
    if (urbanTF) {
      SEL.urban <-
        oco2.lon >= (target.region[5] - dlon) &
        oco2.lon <= (target.region[5] + dlon) &
        oco2.lat >= (target.region[6] - dlat) &
        oco2.lat <= (target.region[6] + dlat)
      tot.urban.count <- length(which(SEL.urban))
      qf.urban.count  <- length(which(SEL.urban & oco2.qf == 0))
      if (oco2.ver == 'b8r')  wl.urban.count <- length(which(SEL.urban & oco2.wl == 0))
      if (oco2.ver == 'b7rb') wl.urban.count <- length(which(SEL.urban & oco2.wl < 15))
    } # end urban TF

    # store results if there are soundings over
    if (tot.count > 0) {

      options(scipen = 999)  # not use exponential express, e.g., 1e10...
      oco.hr <- substr(ncvar_get(oco2.dat, 'sounding_id'), 9, 10)

      tmp.timestr <- paste0(timestr[f], oco.hr)
      uni.timestr <- unique(tmp.timestr[SEL])

      if (urbanTF) {
        tmp <- cbind(as.numeric(uni.timestr), as.numeric(tot.count),
          as.numeric(qf.count), as.numeric(wl.count),
          as.numeric(tot.urban.count), as.numeric(qf.urban.count),
          as.numeric(wl.urban.count))
      } else {
        tmp <- cbind(as.numeric(uni.timestr), as.numeric(tot.count),
          as.numeric(qf.count), as.numeric(wl.count))
      } # end if urbanTF

      result <- rbind(result, tmp)
      nc_close(oco2.dat)

    } else {
      nc_close(oco2.dat)
      next
    } # end if tot.count

  }  # end for f

  result <- as.data.frame(result)

  if (urbanTF) {
    colnames(result) <- c('timestr', 'tot.count', 'qf.count', 'wl.count',
      'tot.urban.count', 'qf.urban.count', 'wl.urban.count')
  } else {
    colnames(result) <- c('timestr', 'tot.count', 'qf.count', 'wl.count')
  } # end if urbanTF

  result <- result[order(result$timestr),]

  return(result)
}
